# INTRODUCTION

Qudaich (queries and unique database alignment inferred by clustering homologs) is a software package for aligning sequences. Qudaich generates the pairwise local alignments between a query dataset against a database. The main design purpose of qudaich is to focus on datasets from next generation sequencing. These datasets generally have hundreds of thousand sequences or more, and the input database will likely contain a large number of sequences. Qudaich is flexible and its algorithmic structure imposes no restriction on the absolute limit of the acceptable read length, but the current version of qudaich allows read lengths &lt;2000 bp. Qudaich can be used to align DNA, translated DNA and protein sequences.

Qudaich performs local sequence alignments in two steps:

1.  Identify the candidate database sequences for each query sequence. The candidate database sequence is the database sequence that gives the best alignment or close enough to the best alignment with the corresponding query sequence.
2.  Generate the optimal alignments between the query sequences and their candidate database sequences using Smith-Waterman-Gotoh algorithm.

Qudaich was written by:

Sajia Akhter, Ph.D. in the [Edwards Bioinformatics Lab](http://edwards.sdsu.edu/research/) and the [Computational Science Research Center](http://www.csrc.sdsu.edu/csrc/) at [San Diego State University](http://www.sdsu.edu/)
and Rob Edwards, Ph.D.

# COPYRIGHT

Qudaich is Copyright 2010-2016 Sajia Akhter and Robert Edwards. It is released under the MIT License, please see the [license file](LICENSE)

# INSTALLATION

1.  Clone this repository using the *Clone or Download* button above.
2.  Compile the code using the command `% make all`. This will result in two executables, *qudaich_alignment* and *qudaich_search_db*

*Qudaich is written in C/C++. So it requires gcc - the GNU project C and C++ compiler - version 4.4.1 or later.*

# QUICK START

## Search The Database

Use qudaich search to find the candidate database sequences:

```% ./qudaich_search_db options```

### Options

* -query Name of the query file (Required)
* -ref Name of database file (Required)
* -prog can be one of the three alignment options (Required):
	* n (nucleotide), 
	* p (protein), 
	* trn (translated nucleotide) 
* -top Number of alignments per query sequence (default 1)
* -freqFile Frequency file Name (default freq.txt)
* -heuristic heuristic options (see below)
	* 1 (default)  or
	* 2

* -h Show command line options

## Generate Alignments

To generate the optimal alignments use the alignment command: 

```% ./qudaich_alignment options```

### Options

* -f Options:
	* all: generate alignments for all query sequences  
	* avg (default): generate alignments for those query sequences whose frequency or sum(lcp) >= average of all query sequences
	* an integer value = generate alignments for those query sequences whose frequency or sum(lcp) >= given integer value

* -freqFile: Name of the frequency file (default: freq.txt) This is the output file generated from `./qudaich_search_db`) 
* -output: Name of output file (default: output\_qudaich.txt)
* -match Match weight (default 1)
* -mismatch Mismatch penalty (default -3)
* -gap\_open Gap opening penalty (default -1)
* -gap\_ext Gap extension penalty (default -2)
* -h Show command line options

# Output


The output file of the alignment is a tab delimited file where the colums are: 

1. query id 
2. database id 
3. query sequence length 
4. identity (not percent identity; you can calculate %identity from column 3 and column 4: %identity = identity/query sequence length * 100) 
5. alignment length 
6. number of gaps
7. query start 
8. query end 
9. reference start 
10. reference end 
11. score


# Heuristic Options

There are two different heuristics that we have developed with qudaich, and you can select them using the -heuristic option. The heuristics are:
Heuristic I: Query *q* has the best alignment with database sequence *d* if the suffixes of *d* are the most frequent closest suffixes in all the query groups containing all the suffixes of *q*.
Heuristic II: Query *q* has the best alignment with database sequence *d*, if heuristic I satisfies and sigma(lcp(suffixes of *q*, suffixes of *d*)) is maximal.


The first heuristic relies on the assertion that if *q* has the best alignment with *d*, the number of common prefix matches (of length at least one) between the suffixes of *q* and the corresponding suffixes of *d* will be the maximum among the number of the common prefix match between the suffixes of *q* and the suffixes of any other database sequence. Since all the suffixes are lexicographically sorted in the suffix array, in most of the cases suffixes of *d* (that have some match with suffixes of *q*) will be either *topDB* or *bottomDB* matches of the query groups containing the corresponding suffixes of *q*.


The second heuristic uses a weighted frequency. Instead of only counting the presence of the DB suffix in *topDB* and *bottomDB*, the longest common prefix (lcp) between the query suffix and the DB suffix is also measured. 

For more information about the choices, see [our paper](http://edwards.sdsu.edu/)




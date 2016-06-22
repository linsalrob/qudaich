# Makefile for the src/ directory. This will be run from the main makefile, so you don't need to run it

CC = gcc
CPP = g++
CFLAGS = -c -O3 -Wall -Wno-char-subscripts -Wno-unused-but-set-variable -Wno-maybe-uninitialized -Wno-unused-variable
LFLAGS = -O3
##CFLAGS = -c -pg -O3 -Wall
##LFLAGS = -pg -O3 -Wall

SOURCEDIR = .
EXECDIR = ../bin

all: |$(EXECDIR) search_dna_hypo1 search_dna_hypo2_top search_pro_hypo2 search_trn_hypo1_top search_dna_hypo1_top search_pro_hypo1 search_pro_hypo2_top search_trn_hypo2 search_dna_hypo2 search_pro_hypo1_top search_trn_hypo1 search_trn_hypo2_top search_trnx_hypo1_top align_dna align_dna_top align_pro align_pro_top align_trn align_trn_top align_trnx_top

search_dna_hypo1: search_db_hypo1_dna_dna.o sais.o
	${CPP} ${LFLAGS} search_db_hypo1_dna_dna.o sais.o -o $(EXECDIR)/search_dna_hypo1
search_dna_hypo1_top: search_db_hypo1_dna_dna_topN.o sais.o
	${CPP} ${LFLAGS} search_db_hypo1_dna_dna_topN.o sais.o -o $(EXECDIR)/search_dna_hypo1_top
search_dna_hypo2: search_db_hypo2_dna_dna.o sais_lcp.o
	${CPP} ${LFLAGS} search_db_hypo2_dna_dna.o sais_lcp.o -o $(EXECDIR)/search_dna_hypo2
search_dna_hypo2_top: search_db_hypo2_dna_dna_topN.o sais_lcp.o
	${CPP} ${LFLAGS} search_db_hypo2_dna_dna_topN.o sais_lcp.o -o $(EXECDIR)/search_dna_hypo2_top
search_pro_hypo1: search_db_hypo1_protein_protein.o sais.o
	${CPP} ${LFLAGS} search_db_hypo1_protein_protein.o sais.o -o $(EXECDIR)/search_pro_hypo1                 
search_pro_hypo1_top: search_db_hypo1_protein_protein_topN.o sais.o    
	${CPP} ${LFLAGS} search_db_hypo1_protein_protein_topN.o sais.o -o $(EXECDIR)/search_pro_hypo1_top
search_pro_hypo2: search_db_hypo2_protein_protein.o sais_lcp.o
	${CPP} ${LFLAGS} search_db_hypo2_protein_protein.o sais_lcp.o -o $(EXECDIR)/search_pro_hypo2
search_pro_hypo2_top: search_db_hypo2_protein_protein_topN.o sais_lcp.o
	${CPP} ${LFLAGS} search_db_hypo2_protein_protein_topN.o sais_lcp.o -o $(EXECDIR)/search_pro_hypo2_top
search_trn_hypo1: search_db_hypo1_translated_translated.o sais.o
	${CPP} ${LFLAGS} search_db_hypo1_translated_translated.o sais.o -o $(EXECDIR)/search_trn_hypo1
search_trn_hypo1_top: search_db_hypo1_translated_translated_topN.o sais.o
	${CPP} ${LFLAGS} search_db_hypo1_translated_translated_topN.o sais.o -o $(EXECDIR)/search_trn_hypo1_top
search_trn_hypo2: search_db_hypo2_translated_translated.o sais_lcp.o
	${CPP} ${LFLAGS} search_db_hypo2_translated_translated.o sais_lcp.o -o $(EXECDIR)/search_trn_hypo2
search_trn_hypo2_top: search_db_hypo2_translated_translated_topN.o sais_lcp.o
	${CPP} ${LFLAGS} search_db_hypo2_translated_translated_topN.o sais_lcp.o -o $(EXECDIR)/search_trn_hypo2_top
search_trnx_hypo1_top: search_db_hypo1_translated_protein_topN.o sais.o
	${CPP} ${LFLAGS} search_db_hypo1_translated_protein_topN.o sais.o -o $(EXECDIR)/search_trnx_hypo1_top
align_dna: sw_affin_dna.o
	${CPP} ${LFLAGS} sw_affin_dna.o -o $(EXECDIR)/align_dna
align_dna_top: sw_affin_dna_topN.o
	${CPP} ${LFLAGS} sw_affin_dna_topN.o -o $(EXECDIR)/align_dna_top
align_pro: sw_affin_protein.o
	${CPP} ${LFLAGS} sw_affin_protein.o -o $(EXECDIR)/align_pro
align_pro_top: sw_affin_protein_topN.o
	${CPP} ${LFLAGS} sw_affin_protein_topN.o -o $(EXECDIR)/align_pro_top
align_trn: sw_affin_translate.o
	${CPP} ${LFLAGS} sw_affin_translate.o -o $(EXECDIR)/align_trn
align_trn_top: sw_affin_translate_topN.o
	${CPP} ${LFLAGS} sw_affin_translate_topN.o -o $(EXECDIR)/align_trn_top
align_trnx_top: sw_affin_translate_protein_topN.o
	${CPP} ${LFLAGS} sw_affin_translate_protein_topN.o -o $(EXECDIR)/align_trnx_top

.c.o: |$(EXECDIR)
	${CC} ${CFLAGS} $<
.cpp.o: |$(EXECDIR)
	${CPP} ${CFLAGS} $<

$(EXECDIR):
	@echo Making $(EXECDIR)
	mkdir -p $(EXECDIR)

.PHONY: clean

clean:
	rm -f $(SOURCEDIR)/*.o $(addprefix $(EXECDIR)/, search_dna_hypo1 search_dna_hypo1_top search_dna_hypo2 search_dna_hypo2_top search_pro_hypo1 search_pro_hypo1_top \
	search_pro_hypo2 search_pro_hypo2_top search_trn_hypo1 search_trn_hypo1_top search_trn_hypo2 search_trn_hypo2_top search_trnx_hypo1_top align_dna \
	align_dna_top align_pro align_pro_top align_trn align_trn_top align_trnx_top)



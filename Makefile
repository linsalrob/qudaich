# Make file for qudaich. This should make qudaich and all the prerequsite programs.

CC = gcc
CPP = g++
CFLAGS = -c -O3 -Wall -Wno-char-subscripts -Wno-unused-but-set-variable -Wno-maybe-uninitialized -Wno-unused-variable
LFLAGS = -O3
##CFLAGS = -c -pg -O3 -Wall
##LFLAGS = -pg -O3 -Wall

SOURCEDIR = src
EXECDIR = bin

SUBDIRS = $(SOURCEDIR)

.PHONY: subdirs $(SUBDIRS) clean


all: |$(EXECDIR) subdirs qudaich_search_db qudaich_alignment

qudaich_search_db: 
	${CPP} ${CFLAGS} $(SOURCEDIR)/run_search_db.cpp -o $(SOURCEDIR)/run_search_db.o
	${CPP} ${LFLAGS} $(SOURCEDIR)/run_search_db.o -o qudaich_search_db
qudaich_alignment: 
	${CPP} ${CFLAGS} $(SOURCEDIR)/run_alignment.cpp -o $(SOURCEDIR)/run_alignment.o
	${CPP} ${LFLAGS} $(SOURCEDIR)/run_alignment.o -o qudaich_alignment

.c.o:
	${CC} ${CFLAGS} $<
.cpp.o:
	${CPP} ${CFLAGS} $<

clean:
	rm -f $(SOURCEDIR)/*.o bin/* qudaich_search_db qudaich_alignment

subdirs: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@

$(EXECDIR):
	@echo Making $(EXECDIR)
	mkdir -p $(EXECDIR)



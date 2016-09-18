VERSION=snapshot
# $Id: inla-Makefile,v 1.21 2009/09/06 16:40:13 hrue Exp $

# USER OPTIONS

# the full path to where to install the library, headerfiles and docs.
# if you use a relative path, you have to modify the Makefile in the
# $(PREFIX)/doc/examples yourself.
PREFIX = /usr/local

# select compilers and optimized compiler-options, add -fopenmp or similar to use OpenMP
CC = gcc
CXX = g++
FC = gfortran
FCEXTRAFLAGS = -fno-second-underscore
FLAGS= -std=gnu99 -O3 -march=native -mfpmath=sse -msse2 -funroll-loops -ftracer -fopenmp -pipe
LD = $(CXX)

# The path to the external libraries: metis, taucs, lapack, blas and zlib, often the same as $(PREFIX)
LEXTPREFIX = $(PREFIX)

# The external libraries to link with
EXTLIBS1 = -L$(PREFIX)/lib -lGMRFLib -L$(LEXTPREFIX)/lib 
EXTLIBS2 = -lgsl -ltaucs -lmetis -llapack -lblas -lgslcblas -lamd -lmuparser -lz -lgfortran
EXTLIBS3 = -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -lm

EXTLIBS = $(EXTLIBS1) $(EXTLIBS2) $(EXTLIBS3)

################################################################################
#
# no changes usually needed below this point
#
################################################################################

SHELL   = /bin/sh

# Only try to update the hgstamp if hg is available.
# Should also check if we are in a repository.
# To use existing HG file, use "make HG= target"
# Would be prettier and more robust with a mercurial update hook.
HG:= $(shell which hg)
HGVERSION:= $(shell if [ -x "$(HG)" ]; then $(HG) parents --template 'hgid: {node|short}  date: {date|date}'; elif [ -e hgstamp ]; then cat hgstamp; else echo 'hgid: unknown  date: unknown' ; fi)

INCL    = -Isrc -I$(PREFIX)/include -DGMRFLib_PATH=\"$(PREFIX)/doc/GMRFLib/data\"
LDFLAGS = $(FLAGS) 
FFLAGS  = $(FLAGS) $(FCEXTRAFLAGS) $(INCL)
CFLAGS  = $(FLAGS) $(INCL) -I$(LEXTPREFIX)/include -DHGVERSION="\"${HGVERSION}\""

INLA = inla
OBJ = src/inla.o src/dictionary.o  src/iniparser.o  src/strlib.o src/my-fix.o src/sha1.o src/spde.o src/spde2.o src/spde3.o src/eval.o src/interpol.o src/re.o src/ar.o src/bessel.o src/pc-priors.o src/R-interface.o src/fgn.o
all: $(INLA)

$(INLA) : $(OBJ)
	$(LD) $(LDFLAGS) -o $@ $^ -L$(PREFIX)/lib  $(EXTLIBS)

install: $(INLA)
	 @umask 022;\
	 mode=644;\
	 for dir in  " " bin doc doc/inla doc/inla/examples; do \
	   test -d $(PREFIX)/$$dir || mkdir $(PREFIX)/$$dir; \
	 done; \
	 install -m755 $(INLA) $(PREFIX)/bin/$(INLA);\
	 install -m755 $(INLA) $(PREFIX)/bin/$(INLA)-$(VERSION);\
	 rsync -auv --no-p --no-o --no-g --chmod=ugo=rwX --delete doc/* $(PREFIX)/doc/inla; \
	 rsync -auv --no-p --no-o --no-g --chmod=ugo=rwX --delete examples/* $(PREFIX)/doc/inla/examples; 

src/TAGS : $(wildcard src/*.c src/*.h)
	etags --output=$@ $^

clean:; -$(RM) $(OBJ) $(INLA)

.PHONY: depend clean uninstall install 





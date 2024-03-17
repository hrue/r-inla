# USER OPTIONS

PREFIX = /usr/local

RLIB_INC = -DINLA_WITH_LIBR -I/usr/include/R -I/usr/share/R/include
RLIB_LIB = -L/usr/lib -lRmath -L/usr/lib/R/lib -lR

# select compilers and optimized compiler-options, add -fopenmp or similar to use OpenMP
CC = gcc
CXX = g++
FC = gfortran
FCEXTRAFLAGS = -fno-second-underscore
FLAGS= -std=gnu99 -O3 -march=native -mfpmath=sse -msse2 -funroll-loops -ftracer -fopenmp -pipe
LD = $(CXX)
INLA_OS = $(shell uname)

# The path to the external libraries: metis, taucs, lapack, blas and zlib, often the same as $(PREFIX)
LEXTPREFIX = $(PREFIX)

# The external libraries to link with
EXTLIBS1 = -L$(PREFIX)/lib -lGMRFLib -L$(LEXTPREFIX)/lib 
EXTLIBS2 = -lgsl -ltaucs -lmetis -llapack -lblas -lgslcblas -lamd -lmuparser -lz -lgfortran -lcrypto
EXTLIBS3 = -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -lm

################################################################################
#
# no changes usually needed below this point
#
################################################################################

SHELL   = /bin/sh

GITCOMMIT:= $(shell git log --pretty=format:"%H - %cd" | head -n1)

EXTLIBS = $(RLIB_LIB) $(EXTLIBS1) $(EXTLIBS2) $(EXTLIBS3)
INCL    = $(RLIB_INC) -Isrc -I$(PREFIX)/include -DGMRFLib_PATH=\"$(PREFIX)/doc/GMRFLib/data\"
LDFLAGS = $(FLAGS) 
FFLAGS  = $(FLAGS) $(FCEXTRAFLAGS) $(INCL)
CFLAGS  = $(FLAGS) $(INCL) -I$(LEXTPREFIX)/include -DGITCOMMIT="\"${GITCOMMIT}\""

INLA = inla
OBJ = src/inla.o src/dictionary.o  src/iniparser.o  src/strlib.o src/my-fix.o src/spde.o \
	src/spde2.o src/spde3.o src/eval.o src/ar.o src/bessel.o src/pc-priors.o \
	src/R-interface.o src/fgn.o src/quantile-regression.o src/my.o src/stochvol.o src/xerbla.o \
	src/tweedie.o src/pc-powerlink.o src/link-gev.o src/cgeneric.o src/libpardiso.o

$(INLA) : $(OBJ)
	$(LD) $(LDFLAGS) -o $@ $(OBJ) -L$(PREFIX)/lib  $(EXTLIBS)

install: $(INLA)
	 @umask 022;\
	 mode=644;\
	 for dir in  " " bin doc doc/inla doc/inla/examples; do \
	   test -d $(PREFIX)/$$dir || mkdir $(PREFIX)/$$dir; \
	 done; \
	 install -m755 $(INLA) $(PREFIX)/bin/$(INLA);\
	 install -m755 $(INLA) $(PREFIX)/bin/$(INLA);\
	 rsync -auv --no-p --no-o --no-g --chmod=ugo=rwX --delete doc/* $(PREFIX)/doc/inla; \
	 rsync -auv --no-p --no-o --no-g --chmod=ugo=rwX --delete examples/* $(PREFIX)/doc/inla/examples; 

clean:; -$(RM) $(OBJ) $(INLA)

.PHONY: depend clean uninstall install 

VERSION=snapshot

# USER OPTIONS

# the full path to where to install the library, headerfiles and docs.
# if you use a relative path, you have to modify the Makefile in the
# $(PREFIX)/doc/examples yourself.
PREFIX = /usr/local

## *** chose either one of these; rgeneric requires the second version ***
ifeq ($(INLA_WITHOUT_LIBR),1)
  RLIB_INC = -DINLA_WITH_RMATH \
		-I/usr/include -I/usr/include/R -I/usr/share/R/include
  RLIB_LIB = -L/usr/lib -lRmath
else
  RLIB_INC = -DINLA_WITH_RMATH -DINLA_LIBR \
		-I/usr/include -I/usr/include/R -I/usr/share/R/include
  RLIB_LIB = -L/usr/lib -lRmath -L/usr/lib/R/lib -lR
endif

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
EXTLIBS2 = -lgsl -ltaucs -lmetis -llapack -lblas -lgslcblas -lamd -lmuparser -lz -lgfortran
EXTLIBS3 = -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -lm

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
HGVERSION:= $(shell if [ -x "$(HG)" -a -d ../.hg ]; then $(HG) parents --template 'hgid: {node|short}  date: {date|date}'; elif [ -e hgstamp ]; then cat hgstamp; else echo 'hgid: unknown  date: unknown' ; fi)

EXTLIBS = $(RLIB_LIB) $(EXTLIBS1) $(EXTLIBS2) $(EXTLIBS3)
INCL    = $(RLIB_INC) -Isrc -I$(PREFIX)/include -DGMRFLib_PATH=\"$(PREFIX)/doc/GMRFLib/data\"
LDFLAGS = $(FLAGS) 
FFLAGS  = $(FLAGS) $(FCEXTRAFLAGS) $(INCL)
CFLAGS  = $(FLAGS) $(INCL) -I$(LEXTPREFIX)/include -DHGVERSION="\"${HGVERSION}\""

#ifeq ($(INLA_OS),Darwin)
#  NOPARDISO_LIB = libpardiso.dylib
#  NOPARDISO_LDFLAGS = -Wall -dynamiclib -undefined suppress -flat_namespace 
#  $(NOPARDISO_LIB) : $(NOPARDISO_SRC)
#	( cd src; $(CC) $(NOPARDISO_LDFLAGS) -o $@ $^ )
#else ifeq ($(INLA_OS),Windows)
#  NOPARDISO_LIB = libpardiso.dll
#  NOPARDISO_CFLAGS = -c
#  NOPARDISO_LDFLAGS = -shared
#  $(NOPARDISO_LIB) : $(NOPARDISO_OBJ)
#	( cd src; $(CC) $(NOPARDISO_LDFLAGS) -o $@ $^ )
#  $(NOPARDISO_OBJ) : $(NOPARDISO_SRC)
#	( cd src; $(CC) $(NOPARDISO_CFLAGS) -c -o $@ $^ )
#else ifeq ($(INLA_OS),Linux)
#  NOPARDISO_LIB = libpardiso.so
#  NOPARDISO_CFLAGS = -fPIC -Wall
#  NOPARDISO_LDFLAGS = -shared
#  $(NOPARDISO_LIB) : $(NOPARDISO_OBJ)
#	( cd src; $(CC) $(NOPARDISO_LDFLAGS) -o $@ $^ )
#  $(NOPARDISO_OBJ) : $(NOPARDISO_SRC)
#	( cd src; $(CC) $(NOPARDISO_CFLAGS) -c -o $@ $^ )
#endif


INLA = inla
OBJ = src/inla.o src/dictionary.o  src/iniparser.o  src/strlib.o src/my-fix.o src/sha1.o src/spde.o src/spde2.o src/spde3.o src/eval.o src/interpol.o src/re.o src/ar.o src/bessel.o src/pc-priors.o src/R-interface.o src/fgn.o src/quantile-regression.o src/my.o src/stochvol.o src/libpardiso.o

$(INLA) : $(OBJ)
	$(LD) $(LDFLAGS) -o $@ $(OBJ) -L$(PREFIX)/lib  $(EXTLIBS)

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

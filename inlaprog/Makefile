# USER OPTIONS

PREFIX = /usr/local

# R's own layout, asked from the R on THIS machine, so one makefile serves
# every device without editing paths per machine. `?=` keeps every existing
# way of overriding: a command-line RLIB_INC=... (what the CI does) or an
# environment variable wins over the derivation. When no R is on the PATH the
# derivation comes back empty and the classic Ubuntu defaults below apply.
R_HOME_DERIVED := $(shell R --vanilla RHOME 2>/dev/null)
ifneq ($(strip $(R_HOME_DERIVED)),)
RLIB_INC ?= -DINLA_WITH_LIBR -I$(shell "$(R_HOME_DERIVED)/bin/Rscript" --vanilla -e 'cat(R.home("include"))' 2>/dev/null)
RLIB_LIB ?=  -L$(R_HOME_DERIVED)/lib -lRmath
else
RLIB_INC ?= -DINLA_WITH_LIBR -I/usr/include/R -I/usr/share/R/include
RLIB_LIB ?= -L/usr/lib -L/usr/lib/R/lib -lRmath
endif

# select compilers and optimized compiler-options, add -fopenmp or similar to use OpenMP
CC := gcc
CXX := g++
FC := gfortran
FCEXTRAFLAGS := -fno-second-underscore
FLAGS := -O3 -march=native -mfpmath=sse -msse2 -funroll-loops -ftracer -fopenmp -pipe
LD := $(CXX)
INLA_OS = $(shell uname)

# The path to the external libraries: metis, lapack, blas and zlib, often the same as $(PREFIX)
LEXTPREFIX = $(PREFIX)

# The external libraries to link with
EXTLIBS1 = -L$(PREFIX)/lib -lGMRFLib -ltaucs -L$(LEXTPREFIX)/lib 
EXTLIBS2 = -lgsl -lmetis -llapack -lblas -lgslcblas -lmuparser -lz -lgfortran -lcrypto
EXTLIBS3 = -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -lm

################################################################################
#
# no changes usually needed below this point
#
################################################################################

SHELL   = /bin/sh

EXTLIBS = $(RLIB_LIB) $(EXTLIBS1) $(EXTLIBS2) $(EXTLIBS3)
INCL    = $(RLIB_INC) -Isrc -I$(PREFIX)/include -DGMRFLib_PATH=\"$(PREFIX)/doc/GMRFLib/data\"
LDFLAGS = $(FLAGS) 
FFLAGS  = $(FLAGS) $(FCEXTRAFLAGS) $(INCL)
CFLAGS  = $(FLAGS) $(INCL) -I$(LEXTPREFIX)/include

INLA = inla
OBJ = src/inla.o src/dictionary.o  src/iniparser.o  src/strlib.o src/my-fix.o src/spde.o \
	src/spde2.o src/spde3.o src/eval.o src/ar.o src/bessel.o src/pc-priors.o \
	src/R-interface.o src/fgn.o src/quantile-regression.o src/my.o src/stochvol.o src/xerbla.o \
	src/tweedie.o src/pc-powerlink.o src/link-gevit.o src/cgeneric.o src/libpardiso.o src/prw2.o \
	src/cgeneric-mapper.o

HEADERS = src/ar.h src/egp-c012.h src/fgn-prior-tables.h src/iniparser.h src/my.h  \
	src/prw2.h src/spde2.h    src/strlib.h src/bessel.h src/eval.h src/fgn-tables-3.h \
	src/inla.h src/param-constr.h src/quantile-regression.h src/spde3.h  src/tweedie.h \
	src/cgeneric.h src/fgn-code.h src/fgn-tables-4.h src/link-gevit.h src/pc-powerlink.h \
	src/R-interface.h src/spde.h src/dictionary.h src/fgn.h src/fgn-tables.h src/my-fix.h \
	src/pc-priors.h src/rmath.h src/stochvol.h

$(INLA) : $(OBJ)
	$(LD) $(LDFLAGS) -o $@ $(OBJ) -L$(PREFIX)/lib  $(EXTLIBS)

src/inla.o: $(HEADERS) \
	src/inla-sys.c src/inla-map-and-link.c src/inla-testit.c src/inla-graph.c src/inla-Qfunc.c \
	src/inla-likelihood.c src/inla-priors.c src/inla-special-functions.c src/inla-output.c \
	src/inla-modes.c src/inla-classic.c src/inla-read.c src/inla-parse.c src/param-constr.c

src/ar.o: $(HEADERS) src/ar.c

src/bessel.o: $(HEADERS) src/bessel.c

src/cgeneric.o: $(HEADERS) src/cgeneric.c

src/cgeneric-demo.o: $(HEADERS) src/cgeneric-demo.c

src/cgeneric-mapper.o: $(HEADERS) src/cgeneric-mapper.c

src/cloglike-demo-beta.o: $(HEADERS) src/cloglike-demo-beta.c

src/cloglike-demo-gaussian.o: $(HEADERS) src/cloglike-demo-gaussian.c

src/cloglike-demo-poisson.o: $(HEADERS) src/cloglike-demo-poisson.c

src/dictionary.o: $(HEADERS) src/dictionary.c

src/eval.o: $(HEADERS) src/eval.c

src/fgn.o: $(HEADERS) src/fgn.c

src/iniparser.o: $(HEADERS) src/iniparser.c

src/libpardiso.o: $(HEADERS) src/libpardiso.c

src/link-gevit.o: $(HEADERS) src/link-gevit.c

src/my.o: $(HEADERS) src/my.c

src/my-fix.o: $(HEADERS) src/my-fix.c

src/param-constr.o: $(HEADERS) src/param-constr.c

src/pc-powerlink.o: $(HEADERS) src/pc-powerlink.c

src/pc-priors.o: $(HEADERS) src/pc-priors.c

src/prw2.o: $(HEADERS) src/prw2.c

src/qinv.o: $(HEADERS) src/qinv.c

src/quantile-regression.o: $(HEADERS) src/quantile-regression.c

src/R-interface.o: $(HEADERS) src/R-interface.c

src/spde2.o: $(HEADERS) src/spde2.c

src/spde3.o: $(HEADERS) src/spde3.c

src/spde.o: $(HEADERS) src/spde.c

src/stochvol.o: $(HEADERS) src/stochvol.c

src/strlib.o: $(HEADERS) src/strlib.c

src/tweedie.o: $(HEADERS) src/tweedie.c

src/view-binary.o: $(HEADERS) src/view-binary.c

src/xerbla.o: $(HEADERS) src/xerbla.c

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

clean:; -$(RM) $(OBJ) $(INLA) $(DEPS)

.PHONY: depend clean uninstall install 

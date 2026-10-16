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
	src/tweedie.o src/tweedie2.o src/pc-powerlink.o src/link-gevit.o src/cgeneric.o src/prw2.o \
	src/cgeneric-mapper.o
OBJFM = src/fast-math/special-functions.o

HEADERS = src/ar.h src/egp-c012.h src/fgn-prior-tables.h src/iniparser.h src/my.h  \
	src/prw2.h src/spde2.h src/strlib.h src/bessel.h src/eval.h src/fgn-tables-3.h \
	src/inla.h src/param-constr.h src/quantile-regression.h src/spde3.h  src/tweedie.h \
	src/cgeneric.h src/fgn-code.h src/fgn-tables-4.h src/link-gevit.h src/pc-powerlink.h \
	src/R-interface.h src/spde.h src/dictionary.h src/fgn.h src/fgn-tables.h src/my-fix.h \
	src/pc-priors.h src/rmath.h src/stochvol.h src/fast-math/special-functions.h \
	src/tweedie2.h src/cgeneric-mapper.h

## -iquote for the same reason as gmrflib/Makefile: ../gmrflib holds an io.h
## that would shadow the system <io.h> on mingw.
FMFLAGS := -DINLA_WITH_CLONE_TARGETS -iquote src -iquote ../gmrflib -fopenmp-simd -O3 -ffast-math
ifneq ($(filter -DINLA_WITH_MKL,$(FFLAGS)),)
FMFLAGS := -DINLA_WITH_MKL $(FMFLAGS)
endif
ifneq ($(filter -DINLA_WITH_FRAMEWORK_ACCELERATE,$(FFLAGS)),)
FMFLAGS := -DINLA_WITH_FRAMEWORK_ACCELERATE $(FMFLAGS)
endif
ifneq ($(filter -DINLA_WITH_ARMPL,$(FFLAGS)),)
FMFLAGS := -DINLA_WITH_ARMPL $(FMFLAGS)
endif
ifneq ($(filter -DINLA_WITH_OPENBLAS,$(FFLAGS)),)
FMFLAGS := -DINLA_WITH_OPENBLAS $(FMFLAGS)
endif
## Same reason as the identical line in gmrflib/Makefile: the define for a
## backend and the -I for its headers arrive together, and copying only the
## define leaves these files unable to find them. special-functions.c includes
## GMRFLib.h, which includes idxval.h, which includes armpl_sparse.h.
FMFLAGS := $(FMFLAGS) $(filter -I%,$(FFLAGS))
src/fast-math/%.o: src/fast-math/%.c
	$(CC) -c $(FMFLAGS) -o $@ $<

src/%.o: src/%.c
	$(CC) -c $(FLAGS) -o $@ $<

$(INLA) : $(OBJ) $(OBJFM)
	$(LD) $(LDFLAGS) -o $@ $(OBJ) $(OBJFM) -L$(PREFIX)/lib $(EXTLIBS)

src/inla.o: $(HEADERS) \
	src/inla-sys.c src/inla-map-and-link.c src/inla-testit.c src/inla-graph.c src/inla-Qfunc.c \
	src/inla-likelihood.c src/inla-priors.c src/fast-math/special-functions.c src/inla-output.c \
	src/inla-modes.c src/inla-classic.c src/inla-read.c src/inla-parse.c src/param-constr.c

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

clean:; -$(RM) $(OBJ) $(OBJFM) $(INLA) $(DEPS)

.PHONY: depend clean uninstall install 

# $Id: GMRFLib-Makefile,v 1.63 2010/01/12 08:55:11 hrue Exp $

# USER OPTIONS

# the full path to where to install the library, headerfiles and docs.
# if you use a relative path, you have to modify the Makefile in the
# $(PREFIX)/doc/examples yourself.
PREFIX = /usr/local

# what to call the library
GMRFLibNAME = GMRFLib

# select compilers and optimized compiler-options. 
CC = gcc
FC = gfortran
FCEXTRAFLAGS = -fno-second-underscore
FLAGS= -std=gnu99 -Wall -O3 -march=native -mfpmath=sse -msse4 -funroll-loops -ftracer -fopenmp -pipe

# The path to the external libraries: metis, taucs, lapack, blas and zlib
LEXTPREFIX = /usr/local
LEXTLIBS = -L$(LEXTPREFIX)/lib 
IEXTLIBS = -I$(LEXTPREFIX)/include -I/usr/include/suitesparse

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

INCL    = -I.. -I$(PREFIX)/include -DGMRFLib_PATH=\"$(PREFIX)/doc/GMRFLib/data\"
LDFLAGS = $(FLAGS) 
FFLAGS  = $(INCL) $(FLAGS) $(FCEXTRAFLAGS)
CFLAGS  = $(INCL) $(FLAGS) $(IEXTLIBS) -DHGVERSION="\"${HGVERSION}\""
LD      = $(FC)
INSTALL = install
LIBNAME = lib$(GMRFLibNAME).a
LIBNAMEG = lib$(GMRFLibNAME)-geo.a

RM    = rm -f
RMR   = rm -rf
RMDIR = rmdir
SED   = sed

LIBOBJ = problem-setup.o lapack-interface.o graph.o error-handler.o acm582.o \
	GMRFLib-fortran.o optimize.o blockupdate.o gdens.o hidden-approx.o \
	distributions.o globals.o wa.o random.o timer.o hash.o density.o \
	smtp-band.o smtp-taucs.o sparse-interface.o rw.o \
	bitmap.o tabulate-Qfunc.o sphere.o io.o approx-inference.o ghq.o \
	utils.o experimental.o graph-edit.o domin-interface.o auxvar.o \
	design.o version.o integrator.o openmp.o hgmrfm.o seasonal.o matern.o \
	bfgs3.o fmesher-io.o smtp-pardiso.o
LIBOBJG = geo.o
HEADERS = blockupdate.h GMRFLib.h  hidden-approx.h optimize.h hash.h \
	distributions.h gdens.h GMRFLibP.h lapack-interface.h timer.h \
	problem-setup.h error-handler.h globals.h graph.h wa.h random.h \
	smtp-band.h smtp-taucs.h sparse-interface.h rw.h \
	bitmap.h hashP.h compatibility.h taucs.h taucs_private.h ghq.h \
	tabulate-Qfunc.h geo.h geo-coefs2.h geo-coefs3.h sphere.h io.h \
	approx-inference.h density.h utils.h experimental.h graph-edit.h \
	domin-interface.h auxvar.h design.h version.h integrator.h openmp.h \
	init.h hgmrfm.h seasonal.h matern.h bfgs3.h fmesher-io.h smtp-pardiso.h
EXAMPLES = examples/Makefile examples/Makefile.in \
	examples/example-blockupdate.c examples/example-graph1.c \
	examples/example-graph2.c examples/example-sample.c \
	examples/example-wa.c examples/germany.data examples/germany.graph \
	examples/graph1.dat examples/wa-example.dat examples/example-rw.c \
	examples/example-geo.c examples/example-qinv.c \
	examples/example-sphere.c examples/example-rng.c \
	examples/example-error-handler.c examples/example-graph-edit.c \
	examples/example-bitmap.c examples/example-auxvar.c \
	examples/example-approx-1.c examples/example-approx-2.c \
	examples/example-approx-3.c examples/poisson.data \
	examples/example-hgmrfm-1.c examples/data_small.dat \
	examples/example-hgmrfm-2.c
MAKEFILESIN = examples/Makefile.in tutorial/Makefile.in

ALL = $(LIBNAME) $(LIBNAMEG) $(MAKEFILESIN)

all : $(ALL)

$(LIBNAME): $(LIBNAME)($(LIBOBJ))
	ranlib $@

$(LIBNAMEG): $(LIBNAMEG)($(LIBOBJG))
	ranlib $@

# !!! Do not change! /FL
DEPDIR = .deps/.dummy
$(DEPDIR):
	-mkdir $(dir $(DEPDIR))
	touch $(DEPDIR)

%.o: %.c $(DEPDIR)
	$(COMPILE.c) -MMD -MF $(dir $(DEPDIR))$*.d.tmp -o $@ $<
	-@$(SED) -e ':more;/\\/N;s/[ ]*\\\r*\n[ ]*/ /;tmore' \
	  -e 's,\($*\.o\)[ :]*\(.*\),$@ : $$\(wildcard \2\),g' \
	  < $(dir $(DEPDIR))$*.d.tmp > $(dir $(DEPDIR))$*.d
	-$(RM) $(dir $(DEPDIR))$*.d.tmp

%.o: %.F $(DEPDIR)
	$(COMPILE.F) -o $@ $<
	-$(RM) $(dir $(DEPDIR))$*.d
	-@grep -E '^#include' < $< \
	  | $(SED) 's!#include "\([^"]*\)"!$@ : $$\(wildcard \1\)!' \
	  > $(dir $(DEPDIR))$*.d

examples/Makefile.in : $(LIBNAME)
	 ( \
	   echo SHELL=$(SHELL); \
	   echo PREFIX=$(PREFIX); \
	   echo GMRFLibNAME=$(GMRFLibNAME);  \
	   echo CC=$(CC); \
	   echo FC=$(FC); \
	   echo FLAGS=$(FLAGS); \
	   echo LEXTLIBS=$(LEXTLIBS); \
	   echo IEXTLIBS=$(IEXTLIBS); \
	 ) > $@

tutorial/Makefile.in : $(LIBNAME)
	 ( \
	   echo SHELL=$(SHELL); \
	   echo PREFIX=$(PREFIX); \
	   echo GMRFLibNAME=$(GMRFLibNAME);  \
	   echo CC=$(CC); \
	   echo FC=$(FC); \
	   echo FLAGS=$(FLAGS); \
	   echo LEXTLIBS=$(LEXTLIBS); \
	   echo IEXTLIBS=$(IEXTLIBS); \
	 ) > $@

doc:
	(cd examples ; make doc )
	$(SED) 's/^\(PROJECT_NUMBER[^=]*=\)\(.*\)/\1 $(shell date)/' \
	  < GMRFLib-Doxyfile \
	  | doxygen -
	mkdir -p doc/html/figs
	mkdir -p doc/tutorial
	cp -fv figs/*.gif doc/html/figs
#	( cd gmrf_tutorial; doxygen )
#	rsync -auv --delete gmrf_tutorial/doxydoc/ doc/tutorial


install: $(ALL)
	 @umask 022;\
	 mode=644;\
	 for dir in  " " lib include include/GMRFLib doc doc/GMRFLib doc/GMRFLib/examples doc/GMRFLib/tutorial doc/GMRFLib/data; do \
	   test -d $(PREFIX)/$$dir || mkdir $(PREFIX)/$$dir; \
	 done; \
	 $(INSTALL) -m$$mode $(LIBNAME) $(PREFIX)/lib; \
	 $(INSTALL) -m$$mode $(LIBNAMEG) $(PREFIX)/lib; \
	 rsync -auv --no-p --no-o --no-g --chmod=ugo=rwX --delete doc/* $(PREFIX)/doc/GMRFLib; \
	 rsync -auv --no-p --no-o --no-g --chmod=ugo=rwX --delete sphere/* $(PREFIX)/doc/GMRFLib/data; \
	 for header in $(HEADERS); do \
	    $(INSTALL) -m$$mode $$header $(PREFIX)/include/GMRFLib; \
	 done; \
	 for ex in $(EXAMPLES); do \
	    $(INSTALL) -m$$mode $$ex $(PREFIX)/doc/GMRFLib/examples; \
	 done;\
	 for ex in tutorial/*; do\
	     $(INSTALL) -m$$mode $$ex $(PREFIX)/doc/GMRFLib/tutorial; \
	 done;\
	 echo; \
	 echo "***	$(LIBNAME) installed     in $(PREFIX)/lib"; \
	 echo "***	$(LIBNAMEG) installed in $(PREFIX)/lib"; \
	 echo "***	documentation is installed in $(PREFIX)/doc/GMRFLib"; \
	 echo "***	examples are installed     in $(PREFIX)/doc/GMRFLib/examples"; \
	 echo "***	tutorial is  installed     in $(PREFIX)/doc/GMRFLib/tutorial"; \
	 echo "***	header files are installed in $(PREFIX)/include/GMRFLib"; \
	 echo ""
uninstall:
	-$(RM) $(PREFIX)/lib/$(LIBNAME)
	-$(RM) $(PREFIX)/lib/$(LIBNAMEG)
	-$(RMR) $(PREFIX)/doc/GMRFLib
	-for header in $(HEADERS);  do $(RM) $(PREFIX)/include/GMRFLib/$$header; done
	-for ex     in $(EXAMPLES); do $(RM) $(PREFIX)/doc/GMRFLib/$$ex; done
	-$(RMDIR) $(PREFIX)/include/GMRFLib

clean:; -$(RM) $(LIBOBJ) $(LIBOBJG) $(LIBNAME) $(LIBNAMEG) $(MAKEFILESIN)
clean-doc:
	-$(RMR) doc/
	-$(RM) example-doxygen-*.txt
clean-deps:
	-$(RMR) $(dir $(DEPDIR))

depend: $(dir $(DEPDIR))dependencies.d
$(dir $(DEPDIR))dependencies.d: Makefile $(DEPDIR) $(wildcard *.c)
	-$(RM) $(dir $(DEPDIR))dependencies.d
	touch $(dir $(DEPDIR))dependencies.d
	for file in *.c *.h ; do \
          if [ -n "`grep HGVERSION $$file`" ]; then \
	    echo $$file | \
	      $(SED) 's/\(.*\)\(\.[hc]\)$$/\1\.o: hgstamp/' \
	      >> $(dir $(DEPDIR))dependencies.d ; \
	  fi ; \
          if [ -n "`grep RCSId.all $$file`" ]; then \
	    echo $$file | \
	      $(SED) 's/\(.*\)\(\.[hc]\)$$/\1\.o: RCSId.all/' \
	      >> $(dir $(DEPDIR))dependencies.d ; \
	  fi \
        done

RCSId.all: hgstamp
	$(SED) 's/\(.*\)/"\1",/' < hgstamp > RCSId.all

hgstamp: dummytarget
	[ -f $@ ] || touch $@
	echo '$(HGVERSION)' | cmp -s $@ - || echo '$(HGVERSION)' > $@
dummytarget: ;

TAGS : $(wildcard *.c *.h)
	etags $^

.PHONY: depend clean clean-deps uninstall install doc dummytarget 

include $(dir $(DEPDIR))dependencies.d $(wildcard $(dir $(DEPDIR))*.d)

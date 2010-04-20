VERSION=snapshot
# $Id: inla-Makefile,v 1.21 2009/09/06 16:40:13 hrue Exp $

# USER OPTIONS

# The path prefix to the external libraries:
# metis, taucs, lapack, blas and zlib
LEXTPREFIX = /usr/local

# External locations of source packages:
URL_METIS = http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-4.0.tar.gz
URL_TAUCS = http://www.tau.ac.il/~stoledo/taucs/2.2/taucs.tgz

# Local filenames to be used for downloaded source packages:
FILE_METIS = metis-4.0.tgz
FILE_TAUCS = taucs-2.2.tgz

# select compilers and optimized compiler-options, add -fopenmp or similar to use OpenMP
CC = gcc
FC = gfortran
FCEXTRAFLAGS = -fno-second-underscore
FLAGS= -std=gnu99 -O3 -march=native -mfpmath=sse -msse2 -funroll-loops -ftracer -fopenmp -pipe
LD = $(CC)


################################################################################
#
# no changes usually needed below this point
#
################################################################################

RM    = rm
RMR   = rm r-f
RMDIR = rmdir
MKDIR = mkdir
MV    = mv
SED   = sed
WGET  = wget

SHELL   = /bin/sh
TMPFILE:= /tmp/$(shell whoami)-inla-extlibs-build.tmp

# LDFLAGS = $(FLAGS) 
# FFLAGS  = $(FLAGS) $(FCEXTRAFLAGS) $(INCL)
# CFLAGS  = $(FLAGS) $(INCL) -I$(LEXTPREFIX)/include

# INLA = inla
# OBJ = src/inla.o src/dictionary.o  src/iniparser.o  src/strlib.o src/my-fix.o src/sha1.o src/sphere.o

METIS = metis-build/metis-4.0/libmetis.a

EXTLIBS_TARGETS = T-metis T-taucs
EXTLIBS_TARGETS_INACTIVE = T-lapack T-atlas T-gsl T-zlib
EXTLIBS_GET = $(patsubst T-%,get-%,$(EXTLIBS_TARGETS))
EXTLIBS_PATCH = $(patsubst T-%,patch-%,$(EXTLIBS_TARGETS))
EXTLIBS_BUILD = $(patsubst T-%,build-%,$(EXTLIBS_TARGETS))
EXTLIBS_INSTALL = $(patsubst T-%,install-%,$(EXTLIBS_TARGETS))
EXTLIBS_GET_ALL = $(EXTLIBS_GET) $(patsubst T-%,get-%,$(EXTLIBS_TARGETS_INCATIVE))
EXTLIBS_PATCH_ALL = $(EXTLIBS_PATCH) $(patsubst T-%,patch-%,$(EXTLIBS_TARGETS_INACTIVE))
EXTLIBS_BUILD_ALL = $(EXTLIBS_BUILD) $(patsubst T-%,build-%,$(EXTLIBS_TARGETS_INACTIVE))
EXTLIBS_INSTALL_ALL = $(EXTLIBS_INSTALL) $(patsubst T-%,install-%,$(EXTLIBS_TARGETS_INACTIVE))

all: build

get: $(EXTLIBS_GET)
unpack: $(EXTLIBS_PATCH)
build: $(EXTLIBS_BUILD)
install: $(EXTLIBS_INSTALL)
get-all: $(EXTLIBS_GET_ALL)
unpack-all: $(EXTLIBS_PATCH_ALL)
build-all: $(EXTLIBS_BUILD_ALL)
install-all: $(EXTLIBS_INSTALL_ALL)

get-metis: $(FILE_METIS)
get-taucs: $(FILE_TAUCS)
$(FILE_METIS):
	if [ ! -f $(FILE_METIS) ]; then \
		$(WGET) -O - $(URL_METIS) > $(FILE_METIS); \
	fi
$(FILE_TAUCS):
	if [ ! -f $(FILE_TAUCS) ]; then \
		$(WGET) -O - $(URL_TAUCS) > $(FILE_TAUCS); \
	fi

patch-metis: $(FILE_METIS)
	$(RMR) metis-build ; \
	mkdir metis-build ; \
	cd metis-build ; \
	tar -zxf ../$(FILE_METIS) ; \
	cd metis-4.0 ; \
	for file in Lib/proto.h Lib/util.c Lib/kmetis.c Lib/kvmetis.c Lib/mkmetis.c ; do \
	  $(SED) 's/log2/log2i/' < $$file | \
	  $(SED) 's/__log2i/metis_log2i/' > $(TMPFILE) ; \
	  $(MV) $(TMPFILE) $$file ; \
	done; \
	cd ../..
	rsync -av --exclude '*~' metis-patch/ metis-build/metis-4.0/
	touch .patch-metis
.patch-metis:
	make patch-metis
build-metis: .patch-metis
	cd metis-build/metis-4.0 ; \
	make ; \
	cd ../..

patch-taucs: $(FILE_TAUCS)
	$(RMR) taucs-build ; \
	mkdir taucs-build ; \
	cd taucs-build ; \
	tar -zxf ../$(FILE_TAUCS) ; \
	patch -p1 -i ../taucs-bugfix.patch ; \
	./configure ; \
	echo "install:\n\tinstall -m755 \$$(DIRLIB)libtaucs\$$(LIBEXT) $(LEXTPREFIX)/lib/libtaucs.a\n" >> makefile ; \
	cd .. ; \
#	rsync -av --exclude '*~' taucs-patch/ taucs-build/
	touch .patch-taucs
.patch-taucs:
	make patch-taucs
build-taucs: .patch-taucs
	cd taucs-build ; \
	make ; \
	cd ../..

$(METIS) :
	make build-metis
install-metis: $(METIS)
	@umask 022;\
	mode=644;\
	for dir in lib; do \
	  test -d $(LEXTPREFIX)/$$dir || mkdir -p $(LEXTPREFIX)/$$dir; \
	done; \
	install -m755 $(METIS) $(LEXTPREFIX)/lib/libmetis.a

install-taucs: build-taucs
	@umask 022;\
	mode=644;\
	for dir in lib; do \
	  test -d $(LEXTPREFIX)/$$dir || mkdir -p $(LEXTPREFIX)/$$dir; \
	done; \
	cd taucs-build ; \
	make install ; \
	cd ..

clean:
	-$(RMR) metis-build
	-$(RMR) taucs-build
	-$(RM) .patch-metis
	-$(RM) .patch-taucs

clean-files:
	-$(RM) $(FILE_METIS) $(FILE_TAUCS)

.PHONY: get build get-all build-all clean $(EXTLIBS_GET_ALL) $(EXTLIBS_PATCH_ALL) $(EXTLIBS_BUILD_ALL)  $(EXTLIBS_INSTALL_ALL) 

#$(INLA) : $(OBJ)
#	$(LD) $(LDFLAGS) -o $@ $^ -L$(PREFIX)/lib  $(EXTLIBS)
#
#install: $(INLA)
#	 @umask 022;\
#	 mode=644;\
#	 for dir in  " " bin doc doc/inla doc/inla/examples; do \
#	   test -d $(PREFIX)/$$dir || mkdir $(PREFIX)/$$dir; \
#	 done; \
#	 install -m755 $(INLA) $(PREFIX)/bin/$(INLA);\
#	 install -m755 $(INLA) $(PREFIX)/bin/$(INLA)-$(VERSION);\
#	 rsync -auv --no-p --no-o --no-g --chmod=ugo=rwX --delete doc/* $(PREFIX#)/doc/inla; \
#	 rsync -auv --no-p --no-o --no-g --chmod=ugo=rwX --delete examples/* $(P#REFIX)/doc/inla/examples; 
#
#clean:; -$(RM) $(OBJ) $(INLA)
#
#.PHONY: depend clean uninstall install





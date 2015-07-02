SHELL = /bin/bash
FLAGS = -std=gnu99 -march=native -mfpmath=sse -msse2 -O3 -ftracer -pipe
##FLAGS = -std=gnu99 -m32 -mtune=generic -mfpmath=sse -msse2 -O3 -ftracer -pipe
PREFIX = $(PWD)/local
PATH := /usr/local/bin:/opt/local/bin:$(PATH)

all :;
	@(echo "";\
	echo "First time only:";\
	echo "";\
	echo "Edit the Makefile, especially the 'FLAGS' variable, then do";\
	echo "";\
	echo "	$(MAKE) download";\
	echo "	$(MAKE) zlib";\
	echo "	$(MAKE) gsl";\
	echo "	$(MAKE) lapack";\
	echo "	$(MAKE) metis";\
	echo "	$(MAKE) taucs";\
	echo "	$(MAKE) amd";\
	echo "	$(MAKE) muparser";\
	echo "	$(MAKE) GMRFLib";\
	echo "	$(MAKE) inla";\
	echo "";\
	echo "If all goes well, you should have a binary of 'inla' at";\
	echo "";\
	echo "	$(PREFIX)/bin/inla";\
	echo "";\
	echo "If you want to include ATLAS libraries, you have to compile them";\
	echo "manually, and put them into";\
	echo "";\
	echo "	$(PREFIX)/lib";\
	echo "";\
	echo "OOPS: rename the ATLAS library liblapack.a into liblapack-atlas.a,";\
	echo "and link also the lacpack-libary with (after lapack-atlas) as the ";\
	echo "ATLAS one does not contain all the routines needed.";\
	echo "";\
	echo "To later update the 'inla' program, do";\
	echo "";\
	echo "	$(MAKE) update";\
	echo "";)

testmac :;  
	@( if [ -f /usr/include/malloc/malloc.h -a ! -f /usr/include/malloc.h ]; then \
		clear;\
		echo "";\
		echo "A soft-link is required for file 'malloc.h'; do";\
		echo "";\
		echo "    sudo ln -s /usr/include/malloc/malloc.h /usr/include/malloc.h";\
		echo "";\
		exit 1;\
	fi;\
	if sed -v 2>&1 | grep -q 'GNU sed'; then \
		dummy=1;\
	else \
		if type -path gsed 2>&1; then \
			clear;\
			echo "";\
			echo "A soft-link is required so that 'sed' is the GNU sed version; do";\
			echo "";\
			echo "    sudo ln -s $(type -path gsed) /usr/local/bin/sed";\
			echo "";\
		else \
			clear;\
			echo "";\
			echo "We need 'sed' to be the GNU-sed version; do";\
			echo "    1. Install MacPorts from www.macports.org";\
			echo "    2. sudo port install gsed";\
			echo "    3. sudo ln -s /opt/local/bin/gsed /opt/local/bin/sed";\
			echo "";\
			exit 1;\
		fi;\
	fi; )

download :; 
	$(MAKE) testmac && ( \
		test -d r-inla && rm -rf r-inla;\
		hg clone https://bitbucket.org/hrue/r-inla;\
		ln -sf r-inla/* .;\
		mkdir -p local/lib;\
		mkdir -p local/bin;\
		mkdir -p local/include; )

zlib :
	cd ./extlibs/zlib-1.2.5; \
	$(MAKE) $(ARGS) -k clean; \
	export CFLAGS="$(FLAGS)";\
	./configure --prefix=$(PREFIX); \
	$(MAKE) $(ARGS); \
	$(MAKE) $(ARGS) install

gsl :
	cd ./extlibs/gsl-1.14; \
	$(MAKE) $(ARGS) -k clean;\
	export CFLAGS="$(FLAGS)";\
	./configure --prefix=$(PREFIX); \
	$(MAKE) $(ARGS); \
	$(MAKE) $(ARGS) install

inla :
	$(MAKE) -C inlaprog PREFIX=$(PREFIX) FLAGS="$(FLAGS) -fopenmp" $(ARGS) -k clean
	$(MAKE) -C inlaprog PREFIX=$(PREFIX) FLAGS="$(FLAGS) -fopenmp" $(ARGS) \
		EXTLIBS2="-lGMRFLib -lgsl -ltaucs -lmetis $(PREFIX)/lib/liblapack.a $(PREFIX)/lib/libblas.a -lgslcblas $(PREFIX)/lib/libz.a $(PREFIX)/lib/libamd.a $(PREFIX)/lib/libmuparser.a -lgfortran"\
		EXTLIBS3=""
	cp -v -f inlaprog/inla $(PREFIX)/bin

lapack :
	$(MAKE) -C lapack-3.2.1 FLAGS="$(FLAGS)" LOADOPTS="$(FLAGS)" $(ARGS) -k clean
	$(MAKE) -C lapack-3.2.1 FLAGS="$(FLAGS)" LOADOPTS="$(FLAGS)" $(ARGS) -k
	cp -v -f lapack-3.2.1/lapack_LINUX.a $(PREFIX)/lib/liblapack.a
	ar d lapack-3.2.1/blas.a lsame.o
	ar d lapack-3.2.1/blas.a xerbla.o
	ar d lapack-3.2.1/blas.a xerbla_array.o
	cp -v -f lapack-3.2.1/blas_LINUX.a $(PREFIX)/lib/libblas.a

GMRFLib :
	$(MAKE) -C gmrflib PREFIX=$(PREFIX) FLAGS="$(FLAGS) -fopenmp" $(ARGS) -k clean
	$(MAKE) -C gmrflib PREFIX=$(PREFIX) FLAGS="$(FLAGS) -fopenmp" $(ARGS)
	$(MAKE) -C gmrflib PREFIX=$(PREFIX) FLAGS="$(FLAGS) -fopenmp" $(ARGS) install

taucs :
	$(MAKE) -C taucs-2.2--my-fix CFLAGS="$(FLAGS)" FFLAGS="$(FLAGS)" $(ARGS) -k clean
	$(MAKE) -C taucs-2.2--my-fix CFLAGS="$(FLAGS)" FFLAGS="$(FLAGS)" $(ARGS)
	cp -v -f taucs-2.2--my-fix/lib/linux/libtaucs.a $(PREFIX)/lib

metis :
	$(MAKE) -C metis-4.0 OPTFLAGS="$(FLAGS)" $(ARGS) -k clean
	$(MAKE) -C metis-4.0 OPTFLAGS="$(FLAGS)" $(ARGS)
	cp -v -f metis-4.0/libmetis.a $(PREFIX)/lib

amd:
	-$(MAKE) -C amd -k clean
	$(MAKE) -C amd FLAGS="$(FLAGS)" PREFIX=$(PREFIX) CC="$(CC)"
	$(MAKE) -C amd FLAGS="$(FLAGS)" PREFIX=$(PREFIX) install

muparser :
        @ cd muparser_v134;\
        ./configure --prefix=$(PREFIX) --disable-shared;\
	$(MAKE) -k CXXFLAGS="$(FLAGS) -DMUPARSER_DLL"; $(MAKE) install


update :; 
	$(MAKE) download;
	$(MAKE) GMRFLib
	$(MAKE) inla

.PHONY: all zlib gsl inla lapack GMRFLib taucs metis download update testmac muparser amd


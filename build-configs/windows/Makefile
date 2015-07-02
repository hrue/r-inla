SHELL = /bin/bash
FLAGS = -std=gnu99 -m32 -mtune=generic -mfpmath=sse -msse -O3 -ftracer -pipe -DWINDOWS -DOSTYPE_win32
PREFIX = $(PWD)/local
PATH := $(PWD)/bin:$(PATH)
ARGS = 

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
	echo "If all goes well, you should have a (static) binary of 'inla' at";\
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


download :; test -d r-inla && rm -rf r-inla;\
	hg clone https://bitbucket.org/hrue/r-inla;\
	ln -sf r-inla/* .;\
	mkdir -p local/lib;\
	mkdir -p local/bin;\
	mkdir -p local/include;

zlib :
	cd zlib-1.2.3; \
	$(MAKE) $(ARGS) -k clean; \
	export CFLAGS="$(FLAGS)";\
	./configure --prefix=$(PREFIX); \
	$(MAKE) $(ARGS); \
	$(MAKE) $(ARGS) install

gsl :
	cd gsl-1.14; \
	$(MAKE) $(ARGS) -k clean;\
	export CFLAGS="$(FLAGS)";\
	./configure --prefix=$(PREFIX); \
	find . -name \*.Po | while read f; do rm -f $$f; touch $$f; done;\
	find . -name \*.Plo | while read f; do rm -f $$f; touch $$f; done;\
	$(MAKE) $(ARGS);\
	$(MAKE) $(ARGS) install

inla :
	$(MAKE) -C inlaprog PREFIX=$(PREFIX) FLAGS="$(FLAGS) -static -fopenmp" $(ARGS) -k clean
	$(MAKE) -C inlaprog PREFIX=$(PREFIX) FLAGS="$(FLAGS) -static -fopenmp" $(ARGS) EXTLIBS3=""
	cp -v -f inlaprog/inla.exe $(PREFIX)/bin

lapack :
	$(MAKE) -C lapack-3.2.1 FLAGS="$(FLAGS)" LOADOPTS="$(FLAGS)" $(ARGS) -k clean
	$(MAKE) -C lapack-3.2.1 FLAGS="$(FLAGS)" LOADOPTS="$(FLAGS)" $(ARGS) -k
	cp -v -f lapack-3.2.1/lapack.a $(PREFIX)/lib/liblapack.a
	cp -v -f lapack-3.2.1/blas.a $(PREFIX)/lib/libblas.a

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

update :; $(MAKE) download;
	$(MAKE) GMRFLib
	$(MAKE) inla

.PHONY: all zlib gsl inla lapack GMRFLib taucs metis download update muparser amd


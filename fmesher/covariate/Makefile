SHELL=bash

USERNAME = $(shell whoami)
ifeq ($(USERNAME),danielsimpson)
INLAPREFIX=/inla/
EIGENPREFIX=/inla/
endif
ifeq ($(USERNAME),finnkrl)
INLAPREFIX=/home/finnkrl/hg/inla/
EIGENPREFIX=/home/finnkrl/kits/eigen/eigen/
endif

covariates: $(wildcard *.cc) $(wildcard *.hh) $(wildcard *.tcc) Makefile
	g++ -O3 -mfpmath=sse -msse2 -funroll-loops -ftracer -pipe \
	-DMAC -DFMESHER_NO_X  -I$(INLAPREFIX)fmesher/src -I$(EIGENPREFIX) \
	-L$(INLAPREFIX)fmesher/src covariates.cc \
	$(INLAPREFIX)fmesher/src/{basis,ioutils,mesh,meshc,predicates,vector,x11utils}.cc \
	-I/inla/extlibs/gsl-1.14/ -Wall -lgsl -lgslcblas \
	-L/inla/extlibs/local/lib/ \
	-o covariates \
	-g

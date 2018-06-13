## Makefile for downloading the source code and building
## a binary package of R-INLA
##
## Run 'make' with no arguments for usage instructions
##
## This Makefile has been tested under Ubuntu 17.10
##
## Requirements:
##   mercurial
##   A full gcc/g++/gfortran installation
##   lib*-dev packages, see target lib-install below
##
## If you already have a cloned inla code repository in a subdirectory 'inla',
## you can skip all 'download-' steps, and instead run
##   make -f inla/build-user/linux/Makefile
##
## Finn Lindgren, 2014-08-24, finn.lindgren@gmail.com
## Finn Lindgren, 2018-01-30, finn.lindgren@gmail.com

SHELL = /bin/bash

## This Makefile is only supported for platform 'linux':
PLATFORM = linux
## Must match the machine architecture (32 or 64):
BITS = 64
## Compiler flags
FLAGS_COMMON = \
	-m64 -Wall -O2 -mtune=generic -mfpmath=sse \
	-msse2 -funroll-loops -ftracer -pipe

## You probably don't need to edit anything below this line

FLAGS_INLA = $(FLAGS_COMMON) -std=gnu99 -fopenmp
FLAGS_FMESHER = $(FLAGS_COMMON) -std=c++11

## For Debian or Ubuntu package installation in lib-install:
APT=aptitude

## Path configuration:
MAINDIR=`pwd`
PREFIX = $(MAINDIR)/local
BINPREFIX = $(PREFIX)/bin/$(PLATFORM)
LIBPREFIX = $(PREFIX)/lib
INCPREFIX = $(PREFIX)/include

## Compiler options:
CC=gcc
CXX=g++
FC=gfortran
LD=$(CXX)
ARGS = CC="$(CC)" CXX="$(CXX)" FC="$(FC)" LD="$(LD)"

all :;
	@(echo "";\
	echo "Edit this Makefile, especially the 'FLAGS_COMMON' variable,";\
	echo "check if you need to run 'sudo $(MAKE) lib-install', then do";\
	echo "";\
	echo "	$(MAKE) lib-install";\
	echo "	$(MAKE) download";\
	echo "	$(MAKE) bin";\
	echo "	$(MAKE) package";\
	echo "";\
	echo "If all goes well, you should have binaries for inla and fmesher at";\
	echo "";\
	echo "	$(BINPREFIX)";\
	echo "";\
	echo "and an INLA_*.tgz binary R package in your working directory.";\
	echo "";\
	echo "To later download and build an updated package, do";\
	echo "";\
	echo "	$(MAKE) update";\
	echo "";\
	echo "To build a package for the JSS 2015 paper, do";\
	echo "";\
	echo "	$(MAKE) download-jss2014 && $(MAKE) bin && $(MAKE) package";\
	echo "";\
	echo "Note: Not all lapack versions contain all the needed routines,";\
	echo "so make sure that /usr/lib/lapack/liblapack.so is selected in";\
	echo "";\
	echo "	sudo update-alternatives --config liblapack.so";\
	echo "";\
	echo "It's possible to link to both lapack and lapack_atlas.";\
	echo "";\
	echo "For x11 enabled fmesher, use lib-install-x11 and bin-x11.";\
	echo "";\
	echo "To deactivate the use of libR, use '$(MAKE) bin INLA_WITHOUT_LIBR=1'.";\
	echo "";)

## Package list
##  updated 2014-08-23
##  updated 2018-05-23: removed optional fmesher libx11-dev (default fmesher builds without x11)
LIBPKG = liblapack-dev libgsl0-dev zlib1g-dev libsuitesparse-dev \
	libmetis-dev libxdmcp-dev libc6-dev libatlas-base-dev \
	libblas-dev libmuparser-dev r-mathlib
LIBPKGX = $(LIBPKG) libx11-dev

lib-install :
	$(APT) install $(LIBPKG)

lib-install-x11 :
	$(APT) install $(LIBPKGX)

init :
	@mkdir -p $(LIBPREFIX)
	@mkdir -p $(BINPREFIX)
	@mkdir -p $(INCPREFIX)

download : init
	@( test -d inla && $(MAKE) download-update ) || \
	( rm -rf inla && hg clone https://bitbucket.org/hrue/r-inla inla )
	@ln -sTf $(PREFIX) inla/local

download-update : init
	@cd inla ; hg pull ; hg update
	@ln -sTf $(PREFIX) inla/local

download-jss2014 :
	@$(MAKE) download
	@cd inla ; hg update -r jss2014
	@ln -sTf $(PREFIX) inla/local

bin :
	$(MAKE) taucs
	$(MAKE) fmesher
	$(MAKE) GMRFLib
	$(MAKE) inla

bin-x11 :
	$(MAKE) taucs
	$(MAKE) fmesher-x11
	$(MAKE) GMRFLib
	$(MAKE) inla

update :
	$(MAKE) download-update
	$(MAKE) bin
	$(MAKE) package


fmesher-x11:
	$(MAKE) -C inla/fmesher PREFIX=$(PREFIX) $(ARGS) -k clean
	$(MAKE) -C inla/fmesher PREFIX=$(PREFIX) FMESHER_USE_X=1 FLAGS="$(FLAGS_FMESHER)" $(ARGS)
	$(MAKE) -C inla/fmesher PREFIX=$(PREFIX) $(ARGS) install
	@rsync -a $(PREFIX)/bin/fmesher $(BINPREFIX)/fmesher$(BITS)
	@rm $(PREFIX)/bin/fmesher

fmesher :
	$(MAKE) -C inla/fmesher PREFIX=$(PREFIX) $(ARGS) -k clean
	$(MAKE) -C inla/fmesher PREFIX=$(PREFIX) FLAGS="$(FLAGS_FMESHER)" $(ARGS)
	$(MAKE) -C inla/fmesher PREFIX=$(PREFIX) $(ARGS) install
	@rsync -a $(PREFIX)/bin/fmesher $(BINPREFIX)/fmesher$(BITS)
	@rm $(PREFIX)/bin/fmesher

inla :
	ln -sTf /usr/include $(INCPREFIX)/muParser
	$(MAKE) -C inla/inlaprog PREFIX=$(PREFIX) FLAGS="$(FLAGS_INLA)" $(ARGS) -k clean
	$(MAKE) -C inla/inlaprog PREFIX=$(PREFIX) FLAGS="$(FLAGS_INLA)" $(ARGS)
	$(MAKE) -C inla/inlaprog PREFIX=$(PREFIX) FLAGS="$(FLAGS_INLA)" $(ARGS) install
	@rsync -a $(PREFIX)/bin/inla $(BINPREFIX)/inla$(BITS)
	@rm $(PREFIX)/bin/inla $(PREFIX)/bin/inla-snapshot

GMRFLib :
	$(MAKE) -C inla/gmrflib PREFIX=$(PREFIX) FLAGS="$(FLAGS_INLA)" $(ARGS) -k clean
	$(MAKE) -C inla/gmrflib PREFIX=$(PREFIX) FLAGS="$(FLAGS_INLA)" $(ARGS)
	$(MAKE) -C inla/gmrflib PREFIX=$(PREFIX) FLAGS="$(FLAGS_INLA)" $(ARGS) install

sync-taucs :
	mkdir -p $(MAINDIR)/tmp
	-rm -rf $(MAINDIR)/tmp/taucs-2.2--my-fix
	cp inla/extlibs/taucs-2.2--my-fix.tgz $(MAINDIR)/tmp/
	cd $(MAINDIR)/tmp ; tar -zxf taucs-2.2--my-fix.tgz

taucs : sync-taucs
	$(MAKE) -C $(MAINDIR)/tmp/taucs-2.2--my-fix CFLAGS="$(FLAGS_INLA)" FFLAGS="$(FLAGS_INLA)" $(ARGS) -k clean
	$(MAKE) -C $(MAINDIR)/tmp/taucs-2.2--my-fix CFLAGS="$(FLAGS_INLA)" FFLAGS="$(FLAGS_INLA)" $(ARGS)
	@cp -v -f $(MAINDIR)/tmp/taucs-2.2--my-fix/lib/linux/libtaucs.a $(PREFIX)/lib

package :
	@rm -f inla/INLA_*.tgz
	cd inla ; utils/build-package-bin $(PREFIX)/bin $(MAINDIR)/tmp

.PHONY: all inla GMRFLib taucs download download-update update fmesher bin init download-jss2014 sync-taucs package fmesher-x11 bin-x11

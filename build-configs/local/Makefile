## Build inla with sTiles on this machine, using the scripts the CI uses.
##
##     make            fetch libstiles, then build -> <repo>/local/bin/inla
##     make inla       build only, libstiles already fetched
##     make stiles     fetch libstiles only
##     make plain      build WITHOUT sTiles
##     make clean      remove the binary and the staged libstiles
##     make distclean  also the objects and the prefix (stale-object fix)
##     make deps       install system packages (sudo, Ubuntu; opt-in)
##
## Installs nothing by default: it assumes the compilers and libraries are
## already present, fetches a prebuilt libstiles from the public releases, and
## builds. `make deps` is separate precisely because it needs sudo.
##
## Every real decision (compiler, CPU target, flags, library versions) stays in
## ci/toolchain.env and ci/build.sh, which the CI reads too. Nothing is
## duplicated here, so this cannot drift from what the released binaries are
## built with. That is what went wrong with build-configs/my-laptop/Makefile,
## which still says -std=gnu99 -msse4 and knows nothing about sTiles or MKL.
##
## Knobs pass straight through, e.g.
##     make JOBS=4
##     make BLAS=armpl
##     make PREFIX=/tmp/inla-test

ROOT  := $(abspath $(dir $(lastword $(MAKEFILE_LIST)))../..)
CI    := $(ROOT)/ci
PREFIX ?= $(ROOT)/local

## Use the compiler this machine has. ci/toolchain.env pins gcc-16 so the
## RELEASED binaries are all built the same way, and refuses to fall back
## quietly; an empty GCC_PREFER is the documented opt-out for a lane that does
## not want the pin (see the comment above GCC_PREFER there). A local build is
## such a lane: requiring gcc-16 on a laptop that has gcc-14 would mean
## installing a compiler before you can compile anything.
##
## The consequence, worth knowing: a binary built here is NOT flag-identical to
## a released one. For timing comparisons against a release, that difference is
## the first thing to suspect.
## GCC_PREFER empty opts out of the pin. GCC_TOOLSET is set unconditionally in
## toolchain.env so it cannot be cleared that way, and the checks refuse when a
## pin is set but no compiler was resolved. Naming CC/CXX/FC satisfies them and
## is honest: a local build should say which compiler it used.
GCC_PREFER ?=
CC  ?= gcc
CXX ?= g++
FC  ?= gfortran

## The released libstiles is built by the pinned gcc, this build uses whatever
## the machine has, and fetch-stiles.sh refuses to mix them. That refusal is
## right for a bundle somebody ships and wrong for a build on your own laptop,
## so opt out here and nowhere else. A release must never set this.
STILES_ALLOW_COMPILER_MISMATCH ?= 1

## No LTO locally. toolchain.env turns it on for releases, and it requires every
## compiler to be the same version: a machine with gcc-14 but gfortran-13, which
## is an ordinary Ubuntu state, fails the link with "bytecode stream generated
## with LTO version 13.1 instead of the expected 14.0". It also makes a local
## rebuild much slower for no benefit you would notice. Pass LTO=1 to opt in.
LTO ?= 0

export PREFIX GCC_PREFER CC CXX FC STILES_ALLOW_COMPILER_MISMATCH LTO

.PHONY: all inla plain stiles deps clean distclean show

all: stiles inla

## libstiles first: ci/build.sh with WITH_STILES=1 expects it staged at
## <prefix>/stiles/{lib,include}, and fails clearly if it is missing.
stiles:
	bash $(CI)/fetch-stiles.sh

inla:
	WITH_STILES=1 bash $(CI)/build.sh

plain:
	WITH_STILES=0 bash $(CI)/build.sh

## Ubuntu only, and it uses sudo, which is why it is not part of `all`.
## Other platforms: ci/deps-fedora.sh, ci/deps-macos.sh, ci/deps-mkl.sh.
deps:
	bash $(CI)/deps-ubuntu.sh

## The binary and the staged libstiles, not the object files: a rebuild after
## this still reuses what the makefiles already compiled.
clean:
	rm -rf $(PREFIX)/bin/inla $(PREFIX)/stiles

## Also the compiled objects and the installed prefix. Needed more often than it
## sounds: make keeps an object that is newer than its source, so one left from
## an older configuration is archived into the new build and the link then fails
## on a symbol nothing provides any more. A taucs object from the PARDISO era
## did exactly that, and the error named PARDISO long after it was removed.
## The objects go through each makefile's own clean, which knows what it built.
## Deleting gmrflib/taucs/obj outright does NOT work: that makefile writes into
## the directory and never creates it, so the next build dies with
## "can't create obj/amdhaf.o". The directory is tracked in git through a
## .dummy file for exactly that reason, so only the objects are removed.
distclean: clean
	-$(MAKE) -C $(ROOT)/gmrflib clean
	-$(MAKE) -C $(ROOT)/inlaprog clean
	rm -f $(ROOT)/gmrflib/taucs/obj/*.o $(ROOT)/gmrflib/taucs/libtaucs.a
	rm -rf $(ROOT)/local $(PREFIX)

show:
	@echo "repo    : $(ROOT)"
	@echo "prefix  : $(PREFIX)"
	@echo "binary  : $(PREFIX)/bin/inla"
	@echo "stiles  : $(PREFIX)/stiles"

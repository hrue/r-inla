SHELL = /bin/bash

default:;
	@echo "TARGETS: all build-configs update-inla update-R update-r-inla.org update-build-yourself www Rinla-local Rinla-install Rinla-testing"

all:;
	$(MAKE) build-configs update-inla update-R update-r-inla.org update-build-yourself

build-configs :;
	./make-builds

update-inla :;
	-tilix --title="laptop-g"  -e "make -C /home/hrue/p/inla/builds/my-laptop-g update"
	-tilix --title="laptop-O2" -e "make -C /home/hrue/p/inla/builds/my-laptop-O2 update"
	-tilix --title="mac64"     -e "ssh -x -t $$INLA_MAC   export INLA_TAG=$$INLA_TAG\; make -C builds/build64 update"
	-tilix --title="$$INLA_LINUX 64"        -e "ssh -p $$INLA_LINUX_PORT -x -t $$INLA_LINUX export INLA_TAG=$$INLA_TAG\; make -C builds/build64 update"
	-tilix --title="$$INLA_LINUX static"    -e "ssh -p $$INLA_LINUX_PORT -x -t $$INLA_LINUX export INLA_TAG=$$INLA_TAG\; make -C builds/build.static update"
	-tilix --title="$$INLA_LINUX windows32" -e "ssh -p $$INLA_LINUX_PORT -x -t $$INLA_LINUX export INLA_TAG=$$INLA_TAG\; make -C builds/windows32 update"
	-tilix --title="$$INLA_LINUX windows64" -e "ssh -p $$INLA_LINUX_PORT -x -t $$INLA_LINUX export INLA_TAG=$$INLA_TAG\; make -C builds/windows64 update"
	@echo -en "\n\nType YES to proceeed: "
	@read ans && test $$ans == "YES" || exit 1
	@echo -e "Proceed with the build...\n\n"

update-R :;
	make -C R updateREPO

update-REPO:;
	make -C R updateREPO

update-r-inla.org :;
	-rsync -auvL --delete r-inla/r-inla.org/. www/r-inla.org/.

update-build-yourself :;
	-rsync -auv --exclude 'mc.tgz' --exclude 'my-laptop.tgz' build-yourself/.  www/build-yourself/.

www :;
	$(MAKE) update-r-inla.org
	-rsync -auv --delete www/. ~/www/www.inla/. 

Rinla-local :;
	make update-inla
	make update-R
	make -C R local

Rinla-local-Ronly :;
	make update-R
	make -C R local

Rinla-install:;
	make update-inla
	make update-R
	make -C R install
	make www
	-do-rsync-www

Rinla-testing :;
	make update-inla
	make update-R
	make -C R testing
	make www
	-do-rsync-www

Rinla-testing-Ronly :;
	make update-REPO
	make -C R testing
	make www
	-do-rsync-www

.PHONY: build-configs update-inla update-R update-REPO update-r-inla.org update-build-yourself all default www Rinla-testing Rinla Rinla-install Rinla-local Rinla-local-Ronly



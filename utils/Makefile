SHELL = /bin/bash 

default:;
	@echo "TARGETS: all build-configs update-inla update-R update-r-inla.org update-build-yourself www Rinla-local Rinla-install Rinla-testing md5"

all:;
	$(MAKE) build-configs update-inla update-R update-r-inla.org update-build-yourself md5

build-configs :;
	./make-builds

update-inla :;
	-tilix --title="mac64"     -e "ssh -x -t $$INLA_MAC         export INLA_TAG=$$INLA_TAG\; make -C builds/build64 update"
	-tilix --title="mac_arm64" -e "ssh -x -t $$INLA_MAC_ARM64   export INLA_TAG=$$INLA_TAG\; make -C builds/arm64   update"
	-tilix --title="$$INLA_LINUX 64"        -e "ssh -p $$INLA_LINUX_PORT -x -t $$INLA_LINUX export INLA_TAG=$$INLA_TAG\; make -C builds/build64 update"
	-tilix --title="$$INLA_WINDOWS windows64" -e "ssh -x -t $$INLA_WINDOWS export INLA_TAG=$$INLA_TAG\; make -C builds/windows64 update"
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

md5 :; ./www/make-md5

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

.PHONY: build-configs update-inla update-R update-REPO update-r-inla.org update-build-yourself all default www Rinla-testing Rinla Rinla-install Rinla-local Rinla-local-Ronly md5



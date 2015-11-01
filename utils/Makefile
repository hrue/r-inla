SHELL = /bin/bash
MAC = hruemac

default:;
	@echo "TARGETS: all build-configs update-inla update-R update-r-inla.org update-build-yourself www Rinla-local Rinla-install Rinla-testing"

all:;
	$(MAKE) build-configs update-inla update-R update-r-inla.org update-build-yourself

build-configs :;
	./make-builds

update-inla :;
	-gnome-terminal --title="LAPTOP"       -e 'nice make -C /home/hrue/p/inla/builds/my-laptop update'
	-gnome-terminal --title="LAPTOP-O3"    -e 'nice make -C /home/hrue/p/inla//builds/my-laptop-O3 update'
	-gnome-terminal --title="LAPTOP-g"     -e 'nice make -C /home/hrue/p/inla//builds/my-laptop-g update'
	-gnome-terminal --title="MAC32"        -e 'ssh -x -t hrue@$(MAC).math.ntnu.no make -C builds/build32 update'
	-gnome-terminal --title="MAC64"        -e 'ssh -x -t hrue@$(MAC).math.ntnu.no make -C builds/build64 update'
	-gnome-terminal --title="MACdebug"     -e 'ssh -x -t hrue@$(MAC).math.ntnu.no make -C builds/debug update'
	-gnome-terminal --title="hrueMAC32"    -e 'ssh -x -t hrue@hruemac.math.ntnu.no make -C builds/build32 update'
	-gnome-terminal --title="hrueMAC64"    -e 'ssh -x -t hrue@hruemac.math.ntnu.no make -C builds/build64 update'
	-gnome-terminal --title="hrueMACdebug" -e 'ssh -x -t hrue@hruemac.math.ntnu.no make -C builds/debug update'
	-gnome-terminal --title="INLA32"       -e 'ssh -x -t hrue@inla.math.ntnu.no make -C builds/build32 update'
	-gnome-terminal --title="INLA64"       -e 'ssh -x -t hrue@inla.math.ntnu.no make -C builds/build64 update'
	-gnome-terminal --title="INLA-O2"      -e 'ssh -x -t hrue@inla.math.ntnu.no make -C builds/build.O2 update'
	-gnome-terminal --title="INLA-O3"      -e 'ssh -x -t hrue@inla.math.ntnu.no make -C builds/build.O3 update'
	-gnome-terminal --title="INLA-g"       -e 'ssh -x -t hrue@inla.math.ntnu.no make -C builds/build.debug update'
	-$(MAKE) -C $(HOME)/p/inla/builds/windows32 update
	-$(MAKE) -C $(HOME)/p/inla/builds/windows64 update
	-gnome-terminal -e 'ssh -x -t hrue@inla.math.ntnu.no sudo install --verbose --mode=755 builds/build.O2/local/bin/inla /usr/local/bin/inla.binary; sudo install --verbose --mode=755 builds/build.debug/local/bin/inla /usr/local/bin/inla.debug.binary; sleep 3'

update-R :;
	make MAC=$(MAC) -C R updateBIN
	make MAC=$(MAC) -C R updateREPO

update-REPO:;
	make MAC=$(MAC) -C R updateREPO

update-r-inla.org :;
	-rsync -auvL --delete google-code/inla/r-inla.org/. www/r-inla.org/.

update-build-yourself :;
	-rsync -auv --exclude 'mc.tgz' --exclude 'my-laptop.tgz' build-yourself/.  www/build-yourself/.

www :;
	$(MAKE) update-r-inla.org
	-rsync -auv --delete www/. ~/www/www.inla/. 

Rinla-local :;
	make update-inla
	make update-R
	make MAC=$(MAC) -C R local

Rinla-local-Ronly :;
	make update-R
	make MAC=$(MAC) -C R local

Rinla-install:;
	make update-inla
	make update-R
	make MAC=$(MAC) -C R install
	make www
	-do-rsync-www
	-gnome-terminal -e 'ssh -x -t hrue@inla.math.ntnu.no "source .bash_profile; source .bashrc; sudo ls; cd; cd tmp; test -f INLA.tgz && rm INLA.tgz; wget http://www.math.ntnu.no/inla/binaries/INLA.tgz; sudo R CMD INSTALL --library=/usr/local/lib/R/site-library INLA.tgz; sudo ~hrue/bin/fix-permissions /usr/local/lib/R/site-library"'

Rinla-testing :;
	make update-inla
	make update-R
	make MAC=$(MAC) -C R testing
	make www
	-do-rsync-www
	-gnome-terminal -e 'ssh -x -t hrue@inla.math.ntnu.no "source .bash_profile; source .bashrc; sudo ls; cd; cd tmp; test -f INLA.tgz && rm INLA.tgz; wget http://www.math.ntnu.no/inla/binaries/testing/INLA.tgz; sudo R CMD INSTALL --library=/usr/local/lib/R/site-library INLA.tgz; sudo ~hrue/bin/fix-permissions /usr/local/lib/R/site-library"'

Rinla-testing-Ronly :;
	make update-REPO
	make MAC=$(MAC) -C R testing
	make www
	-do-rsync-www
	-gnome-terminal -e 'ssh -x -t hrue@inla.math.ntnu.no "source .bash_profile; source .bashrc; sudo ls; cd; cd tmp; test -f INLA.tgz && rm INLA.tgz; wget http://www.math.ntnu.no/inla/binaries/testing/INLA.tgz; sudo R CMD INSTALL --library=/usr/local/lib/R/site-library INLA.tgz; sudo ~hrue/bin/fix-permissions /usr/local/lib/R/site-library"'

.PHONY: build-configs update-inla update-R update-REPO update-r-inla.org update-build-yourself all default www Rinla-testing Rinla Rinla-install Rinla-local Rinla-local-Ronly



SHELL = /bin/bash

local :;
	bash ./install.me 0

install :;
	bash ./install.me 1

testing :;
	bash ./install.me 2

check :; 
	R CMD check --no-install --no-examples INLA

updateBIN:;
	## the output of this command must be copied
	-for f in $$(ssh $(MAC).math.ntnu.no bin/my.otool builds/build32/local/bin/inla); do echo sync $$f; rsync -auv hrue@$(MAC).math.ntnu.no:$$f extern/mac32; done
	-for f in $$(ssh $(MAC).math.ntnu.no bin/my.otool builds/build64/local/bin/inla); do echo sync $$f; rsync -auv hrue@$(MAC).math.ntnu.no:$$f extern/mac64; done
	-for f in $$(ssh $(MAC).math.ntnu.no bin/my.otool builds/build32/local/bin/fmesher); do echo sync $$f; rsync -auv hrue@$(MAC).math.ntnu.no:$$f extern/mac32; done
	-for f in $$(ssh $(MAC).math.ntnu.no bin/my.otool builds/build64/local/bin/fmesher); do echo sync $$f; rsync -auv hrue@$(MAC).math.ntnu.no:$$f extern/mac64; done
	-rsync -auv hrue@$(MAC).math.ntnu.no:builds/build32/local/bin/inla extern/mac32
	-rsync -auv hrue@$(MAC).math.ntnu.no:builds/build32/local/bin/fmesher extern/mac32
	-rsync -auv hrue@$(MAC).math.ntnu.no:builds/build64/local/bin/inla extern/mac64
	-rsync -auv hrue@$(MAC).math.ntnu.no:builds/build64/local/bin/fmesher extern/mac64
	-rsync -auv hrue@inla.math.ntnu.no:builds/build32/local/bin/inla extern/inla32.linux
	-rsync -auv hrue@inla.math.ntnu.no:builds/build32/local/bin/fmesher extern/fmesher32.linux
	-rsync -auv hrue@inla.math.ntnu.no:builds/build64/local/bin/inla extern/inla64.linux
	-rsync -auv hrue@inla.math.ntnu.no:builds/build64/local/bin/fmesher extern/fmesher64.linux

updateREPO:;
	## the output of this command must be copied
	-hg --cwd r-inla pull; hg --cwd r-inla update

.PHONY: updateBIN updateREPO local install testing check

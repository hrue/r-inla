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
	./updateBin

updateREPO:;
	## the output of this command must be copied
	-hg --cwd r-inla pull; hg --cwd r-inla update -r $$INLA_TAG

.PHONY: updateBIN updateREPO local install testing check

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
	cd r-inla && git checkout devel && git branch -d branch-$$INLA_TAG; git pull --all --tags && git checkout $$INLA_TAG -b branch-$$INLA_TAG

.PHONY: updateBIN updateREPO local install testing check

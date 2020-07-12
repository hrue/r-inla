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
	cd r-inla; git reset --hard; git checkout devel; git branch -D tag-$$INLA_TAG||true; git pull --all --prune; git checkout tags/$$INLA_TAG -b tag-$$INLA_TAG

.PHONY: updateBIN updateREPO local install testing check

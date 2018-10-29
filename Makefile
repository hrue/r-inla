## keep track of updates needed to be done

all:; \
	## update the TAGS-file and the man-files that are generated
	## automatically 
	-make -C rinla/R 
	## update the doc extracted from the R-files
	-make -C r-inla.org/doc
	## update the TAGS-file
	-make -C gmrflib TAGS
	## update the TAGS-file
	-make -C inlaprog src/TAGS
	##
	-make doc-links


## The following finds all .tex files in
##   r-inla.org/doc/{prior,latent,likelihood}
## and adds links in rinla/inst/doc/... to the corresponding .pdf files
## In the past it redirected to
##   ../../../../../../google-code/inla/
## Now it redirects to
##   ../../../../
## which keeps it inside the mercurial clone.
doc-links:
	@echo "Building documentation links from r-inla.org to rinla."
	@for dir in prior latent likelihood link; do \
	  find "r-inla.org/doc/$$dir" -name \*.tex | \
		sed "s!r-inla.org/doc/\(.*\)\.tex!ln -sf ../../../../r-inla.org/doc/\1.pdf rinla/inst/doc/\1.pdf!" | sh -e ;\
	  hg status r-inla.org/doc/$$dir/*.tex | cat;\
	done
	hg status rinla/inst/doc/ | cat
	@echo "Documentation link building finished."

## Build a INLA-package without the binaries
INLA-package:
	utils/build.package

## Check the latest built package from INLA-package:
R = R
## R = ~/R-dev/bin/R
#FILENAME := $(shell ls -rt ~/tmp/.BuildINLA/INLA_*.tgz | grep -v _R_ | tail -1 )
INLA-package-check:
	@echo Checking $(FILENAME)
	$(R) --vanilla CMD check --no-examples $(FILENAME)

##
.PHONY: doc-links INLA-package all INLA-package-check

## keep track of updates needed to be done

all:; \
	make tags
	-make -C rinla/R
	-make -C r-inla.org/doc
	-make doc-links


## The following finds all .tex files in
##   r-inla.org/doc/{prior,latent,likelihood}
## and adds links in rinla/inst/documentation/... to the corresponding .pdf files
## In the past it redirected to
##   ../../../../../../google-code/inla/
## Now it redirects to
##   ../../../../
## which keeps it inside the mercurial clone.
doc-links:
	@echo "Building documentation links from r-inla.org to rinla."
	@for dir in prior latent likelihood link; do \
	  find "r-inla.org/doc/$$dir" -name \*.tex | \
	        grep -v "old-stuff" | \
		sed "s!r-inla.org/doc/\(.*\)\.tex!ln -sf ../../../../r-inla.org/doc/\1.pdf rinla/inst/documentation/\1.pdf!" | sh -e ;\
	done; 
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

tags:; gtags; ## htags --line-number=5

##
.PHONY: doc-links INLA-package all INLA-package-check tags

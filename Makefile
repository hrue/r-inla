## keep track of updates needed to be done

all:; \
	## update the doc extracted from the R-files
	-make -C r-inla.org/doc
	## update the TAGS-file and the man-files that are generated
	## automatically 
	-make -C rinla/R 
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
	@for dir in prior latent likelihood; do \
	  find "r-inla.org/doc/$$dir" -name \*.tex | \
		sed "s!r-inla.org/doc/\(.*\)\.tex!ln -sf ../../../../r-inla.org/doc/\1.pdf rinla/inst/doc/\1.pdf!" | sh -e ;\
	  hg status r-inla.org/doc/$$dir/*.tex ;\
	done
	hg status rinla/inst/doc/
	@echo "Documentation link building finished."

## Build a INLA-package without the binaries
INLA-package:
	utils/build.package


##
.PHONY: doc-links INLA-package all

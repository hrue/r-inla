# By default assumes that
#   R
# calls the working R binary.
# Override the R binary with
#   make R=path-to-R-binary
# 
R=R

RR=$(R) --vanilla --no-echo

all: TAGS ../NAMESPACE

.PHONY : all clean models_doc \
	check_Rd_pdf clean_Rd_pdf

namespace: ../NAMESPACE

clean :
	rm -f make-man-pages.timestamp
	rm -f ../NAMESPACE.ours
	rm -f ../NAMESPACE.roxygen
	$(MAKE) clean_Rd_pdf

TAGS : $(wildcard *.R)
	$(RR) -e 'library(utils); rtags(path=".", ofile="TAGS")'


## Replacement for old models.Rd generator. This generates roxygen documentation
## that models.Rd is then generated from.
models_doc: ../R/models_documentation.R
../R/models_documentation.R: inlaEnv.R models.R models-generate.R utils.R version.R
	{ \
	  echo 'options(keep.source=TRUE)';\
	  echo 'source("inlaEnv.R")';\
	  echo 'source("models.R")';\
	  echo 'source("models-generate.R")';\
	  echo 'source("utils.R")';\
	  echo 'source("version.R")';\
	  echo 'inla.models.generate.roxygen("../R/models_documentation.R")';\
	  echo 'styler::style_file("../R/models_documentation.R")';\
	} | $(RR)


## creating roxygen documentation also updates the NAMESPACE file

../NAMESPACE: $(wildcard *.R) ../R/models_documentation.R ../DESCRIPTION Makefile
	-test -f ../NAMESPACE && cp ../NAMESPACE ../NAMESPACE.backup
	echo 'testthat_print <- function(...) print(...) ; devtools::document()' | $(RR)

# Generate the pdf files for the manual one by one to
# get easier information about location of latex issues 
check_Rd_pdf: $(patsubst %.Rd,%.pdf,$(wildcard ../man/*.Rd))

../man/%.pdf: ../man/%.Rd
	$(RR) CMD Rd2pdf --force --no-preview --output=$@ $<

clean_Rd_pdf:
	rm -f ../man/*.pdf

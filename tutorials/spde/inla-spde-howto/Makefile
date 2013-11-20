RnwSOURCES = $(wildcard *.Rnw)
texDepends = $(RnwSOURCES:.Rnw=.tex) 

all : inla-spde-howto.pdf

inla-spde-howto.pdf : $(texDepends)
	@echo generating the PDF file
	latexmk inla-spde-howto.tex -pdf -recorder </dev/null

%.tex : %.Rnw 
	@echo generating file $@
	@R --vanilla -e 'library(knitr); knit("$^")' </dev/null

install:; distill inla-spde-howto.pdf
	cp -v inla-spde-howto.pdf ../../../r-inla.org/tutorials/spde/inla-spde-howto.pdf

clean :; latexmk -C

.PHONY: clean



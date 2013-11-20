RnwSOURCES = $(wildcard *.Rnw)
texDepends = $(RnwSOURCES:.Rnw=.tex) 

all : inla-spde-howto.pdf

inla-spde-howto.pdf : $(texDepends)
	@echo generating the PDF file
	latexmk inla-spde-howto.tex -pdf -recorder </dev/null

%.tex : %.Rnw 
	@echo generating file $@
	@R --vanilla -e 'library(knitr); knit("$^")' </dev/null

clean :; latexmk -C

.PHONY: clean



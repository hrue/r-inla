RnwSOURCES = $(wildcard *.Rnw)
texDepends = updates.tex spde-tutorial.tex spde-tutorial-introduction.tex $(RnwSOURCES:.Rnw=.tex) spde-tutorial.bib

all : spde-tutorial.pdf

spde-tutorial.pdf : $(texDepends)
	@echo generating spde-tutorial.pdf
	latexmk spde-tutorial.tex -pdf -recorder -bibtex </dev/null

%.tex : %.Rnw 
	@echo generating file $@
	@R --vanilla -e 'Sweave("$^")' </dev/null

install:; distill spde-tutorial.pdf
	cp spde-tutorial.pdf ../../../inla/r-inla.org/tutorials/spde/spde-tutorial.pdf

clean :; latexmk -C

.PHONY: clean install



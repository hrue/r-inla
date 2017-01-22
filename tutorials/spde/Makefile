RnwSOURCES = $(wildcard *.Rnw)
texDepends = updates.tex spde-tutorial.tex spde-tutorial-introduction.tex spde-tutorial-spde-results.tex $(RnwSOURCES:.Rnw=.tex) spde-tutorial.bib 

all : pdf html

pdf : $(texDepends)
	@echo generating spde-tutorial.pdf
	latexmk spde-tutorial.tex -pdf -recorder -bibtex </dev/null

html : $(texDepends)
	@echo generating the html output
	htlatex spde-tutorial.tex "html,frames,2"

%.tex : %.Rnw 
	@echo generating file $@
	@R --vanilla -e 'Sweave("$^")' </dev/null

spde-tutorial-spde-introduction.tex : spde-tutorial-functions.R spde-tutorial-spde-introduction.Rnw
	@echo generating file spde-tutorial-spde-introduction.Rnw
	@R --vanilla -e 'Sweave("spde-tutorial-spde-introduction.Rnw")' </dev/null


clean :; latexmk -C

.PHONY: clean



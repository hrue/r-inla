RnwSOURCES = $(wildcard *.Rnw)
texDepends = spde-tutorial.tex spde-tutorial-introduction.tex $(RnwSOURCES:.Rnw=.tex) 

all : spde-tutorial.pdf

spde-tutorial.pdf : $(texDepends)
	@echo generating spde-tutorial.pdf
	latexmk spde-tutorial.tex -pdf -recorder -bibtex </dev/null

%.tex : %.Rnw 
	@echo generating file $@
	@R --vanilla -e 'Sweave("$^")' </dev/null

clean :; latexmk -C

.PHONY: clean



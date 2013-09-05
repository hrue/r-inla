RnwSOURCES = $(wildcard *.Rnw)
texDepends = $(RnwSOURCES:.Rnw=.tex) shs.bib

all : tutorial_igmrfs.tex

tutorial_igmrfs.pdf : $(texDepends)
	@echo generating tutorial_igmrfs.pdf
	latexmk tutorial_igmrfs.tex -pdf -recorder -bibtex </dev/null

%.tex : %.Rnw 
	@echo generating file $@
	@R --vanilla -e 'Sweave("$^")' </dev/null

install :;
	cp -v tutorial_igmrfs.pdf ../../../r-inla.org/tutorials/inla/scale.model/scale-model-tutorial.pdf

clean :; latexmk -C

.PHONY: clean install



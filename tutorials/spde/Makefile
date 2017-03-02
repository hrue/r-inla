RnwSOURCES = $(wildcard *.Rnw)
texDepends = updates.tex spde-tutorial.tex $(RnwSOURCES:.Rnw=.tex) spde-tutorial.bib 

all : pdf html 

pdf : $(texDepends)
	@echo generating spde-tutorial.pdf
	latexmk spde-tutorial.tex -pdf -recorder -bibtex </dev/null

html : $(texDepends)
	@echo generating the html output
	./pdfepsconv.sh
	make4ht -d html spde-tutorial.tex 'html,frames,2'
	mv html/spde-tutorial.html html/index.html
	rm *.html *.png

%.R : %.Rnw
	@R --vanilla -e 'knitr::purl("$^")' </dev/null   
	mv *.R R/

%.tex : %.Rnw 
	@echo generating file $@
	@R --vanilla -e 'library(knitr); knit("$^")' </dev/null

clean :; latexmk -C

.PHONY: clean



RnwSOURCES = $(wildcard *.Rnw)
texDepends = updates.tex spde-tutorial.tex $(RnwSOURCES:.Rnw=.tex) spde-tutorial.bib 
Rdepends = $(RnwSOURCES:.Rnw=.R) 

all : pdf html install

install: spde-tutorial.pdf html/index.html
	@echo coping the files to r-inla.org
	@cp spde-tutorial.pdf ~/bitbucket/r-inla/r-inla.org/tutorials/spde/
	@cp html/* ~/bitbucket/r-inla/r-inla.org/tutorials/spde/html/
	@cp R/* ~/bitbucket/r-inla/r-inla.org/tutorials/spde/R/

pdf : $(texDepends) ##$(Rdepends)
	@echo generating spde-tutorial.pdf
	@latexmk spde-tutorial.tex -pdf -recorder -bibtex </dev/null

html : $(texDepends) ## $(Rdepends)
	@echo generating the html output
	@./pdfepsconv.sh
	@make4ht -d html spde-tutorial.tex 'html,frames,2'
	@mv html/spde-tutorial.html html/index.html
	@rm *.html *.png

%.tex : %.Rnw
	@echo generating file $(@:.tex=.R)
	@R --vanilla -e 'knitr::purl("$^")' </dev/null   
	@mv $(@:.tex=.R) R/
	@echo generating file $@
	@R --vanilla -e 'library(knitr); knit("$^")' </dev/null

clean :; latexmk -C

.PHONY: clean



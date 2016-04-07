RmdSOURCES = $(wildcard *.Rmd)
mdDepends = $(RmdSOURCES:.Rmd=.md) 
htmlDepends = $(mdDepends:.md=.html)
texDepends = SPDE1d.tex SPDEhowto.tex
pdfDepends = SPDE1d.pdf SPDEhowto.pdf

all : igmrfs isbaspde spde-tutorial html pdf

igmrfs : 
	cp ../../tutorials/inla/scale.model/demo.R .
	cp ../../tutorials/inla/scale.model/tutorial_igmrfs.Rnw .
	cp ../../tutorials/inla/scale.model/tutorial_igmrfs.pdf .

spde-tutorial : 
	cp ../../tutorials/spde/spde-tutorial.pdf .

isbaspde : 
	cp ../../tutorials/spde-isba/isbaspde.pdf .

%.md : %.Rmd
	@echo generating file $@
	@R --vanilla -e 'knitr::knit("$^")' </dev/null

%.html : %.md 
	@echo generating file $@
	@Rscript -e "require(markdown); markdownToHTML('$^', '$@', options=c('use_xhtml', 'base64_images'))"

html : $(htmlDepends)

%.pdf : %.md
	@echo generating file $@
	@/usr/lib/rstudio/bin/pandoc/pandoc -V geometry:margin=0.9in -s $^ -o $@ 

tex : $(texDepends)
pdf : $(pdfDepends) 

clean :; latexmk -C

.PHONY: clean


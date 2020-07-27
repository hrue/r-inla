RMDs := $(wildcard *.Rmd)
RNWs := $(wildcard *.Rnw)
PDFs := $(RMDs:.Rmd=.pdf) $(RNWs:.Rnw=.pdf)
HTMLs := $(RMDs:.Rmd=.html) 

all : $(PDFs) $(HTMLs)

%.pdf : %.Rmd
	R --vanilla -e 'library(rmarkdown); render("'$<'", c("pdf_document"))'

%.html : %.Rmd
	R --vanilla -e 'library(rmarkdown); render("'$<'", c("html_document"))'

%.pdf : %.Rnw
	R --vanilla -e 'library(knitr); knit2pdf("'$<'")'


RMDs := $(wildcard *.Rmd)
RNWs := $(wildcard *.Rnw)
PDFs := $(RMDs:.Rmd=.pdf) $(RNWs:.Rnw=.pdf)

all : $(PDFs)

%.pdf : %.Rmd
	R --vanilla -e 'library(rmarkdown); render("'$<'", pdf_document())'

%.pdf : %.Rnw
	R --vanilla -e 'library(knitr); knit2pdf("'$<'")'

compress :;  pdfcompress -v -O $(wildcard *.pdf)

.PHONY : compress

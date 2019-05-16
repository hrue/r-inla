RMDs := $(wildcard *.Rmd)
RNWs := $(wildcard *.Rnw)
PDFs := $(RMDs:.Rmd=.pdf) $(RNWs:.Rnw=.pdf)

all : $(PDFs)

%.pdf : %.Rmd
	R --vanilla -e 'library(rmarkdown); render("'$<'", "all")'

%.pdf : %.Rnw
	R --vanilla -e 'library(knitr); knit2pdf("'$<'")'


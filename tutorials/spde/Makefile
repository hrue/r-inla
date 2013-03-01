all : spde-tutorial.pdf

%.tex : %.Rnw
	R CMD BATCH runCacheSweaveSPDEtutorial.R /dev/stdout </dev/null

%.pdf : %.tex
	latexmk -pdf -recorder -bibtex $< </dev/null

clean :; latexmk -C

.PHONY: clean



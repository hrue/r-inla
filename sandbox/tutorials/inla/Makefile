all : main.pdf

%.pdf : %.tex
	latexmk -pdf -recorder -bibtex $< </dev/null

%.tex : %.Rnw
	R CMD BATCH runCacheSweave.R /dev/stdout </dev/null

clean :; latexmk -C

.PHONY: clean

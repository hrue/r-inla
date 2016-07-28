RmdSources = $(wildcard *.Rmd)
RnwSources = $(wildcard *.Rnw)
PANDOC = /usr/lib/rstudio/bin/pandoc/pandoc
PANDOC = pandoc

all : mdtutorials nwtutorials

%.pdf : %.Rmd
	@echo generating .md file
	@R --vanilla -e 'knitr::knit("$^")' </dev/null
	@echo generating the $@
	@$(PANDOC) -V geometry:margin=0.9in --latex-engine=xelatex -s $(^:.Rmd=.md) -o $@
	@echo generating the HTML
	@$(PANDOC) $(^:.Rmd=.md) --from markdown+autolink_bare_uris+ascii_identifiers+tex_math_single_backslash-implicit_figures -t html -s -o $(@:.pdf=.html) --smart --email-obfuscation none --self-contained --standalone --section-divs ##--mathjax  --variable 'mathjax-url:https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML' 
	@rm figure* -rf 
	@rm *.md

%.tex : %.Rnw 
	@echo generating file $@
	@R --vanilla -e 'Sweave("$^")' </dev/null

%.pdf : %.tex
	@latexmk -pdf $^ </dev/null

mdtutorials : $(RmdSources:.Rmd=.pdf)

nwtutorials : $(RnwSources:.Rnw=.pdf)

clean :; latexmk -C

.PHONY: clean


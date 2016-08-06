TARGETS = master.pdf
DEPS_DIR = .deps
LATEXMK = latexmk -recorder -use-make -deps \
        -e 'warn qq(In Makefile, turn off custom dependencies\n);' \
        -e '@cus_dep_list = ();' \
        -e 'show_cus_dep();'

all : $(TARGETS)

$(foreach file,$(TARGETS),$(eval -include $(DEPS_DIR)/$(file)P))

$(DEPS_DIR) :
	mkdir $@

%.pdf : %.tex
	if [ ! -e $(DEPS_DIR) ]; then mkdir $(DEPS_DIR); fi
	$(LATEXMK) -pdf -dvi- -ps- -deps-out=$(DEPS_DIR)/$@P $<


%.tex : %.Rnw
	R -e 'library(INLA); library(knitr); knit("'$<'")'

clean :; latexmk -CA

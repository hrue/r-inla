Makefile.in :
	( cd .. ; make examples/Makefile.in )
include Makefile.in

INCL    = -I$(PREFIX)/include $(IEXTLIBS)
LDFLAGS = $(FLAGS) $(INCL)
FFLAGS  = $(FLAGS) $(INCL)
CFLAGS  = $(FLAGS) $(INCL)
LD      = $(CC)

PROGS = example-blockupdate \
	example-blockupdate-with-error-control \
	example-graph1 example-graph2 example-sample \
	example-wa example-rw example-geo example-qinv \
	example-sphere example-rng example-error-handler example-graph-edit \
	example-bitmap example-auxvar example-approx-1 example-approx-2 example-approx-3 \
	example-hgmrfm-1 example-hgmrfm-2 example-matern2d
PROGS_MANUAL = $(PROGS) example-blockupdate-2

EXAMPLES_MANUAL = $(PROGS_MANUAL:example-%=../example-doxygen-%.txt)

# add `-lcompat' for FREEBSD

EXTLIBS = $(LEXTLIBS) -lgsl -ltaucs -lmetis -llapack -lblas -lgslcblas -lz -lgfortran -lm

all : $(PROGS) 

example-auxvar : example-auxvar.o
	$(LD) $(FLAGS) -o $@ $< -L$(PREFIX)/lib -l$(GMRFLibNAME) $(EXTLIBS)

example-bitmap : example-bitmap.o
	$(LD) $(FLAGS) -o $@ $< -L$(PREFIX)/lib -l$(GMRFLibNAME) $(EXTLIBS)

example-graph-edit : example-graph-edit.o
	$(LD) $(FLAGS) -o $@ $< -L$(PREFIX)/lib -l$(GMRFLibNAME) $(EXTLIBS)

example-rng : example-rng.o
	$(LD) $(FLAGS) -o $@ $< -L$(PREFIX)/lib -l$(GMRFLibNAME) $(EXTLIBS)

example-blockupdate : example-blockupdate.o
	$(LD) $(FLAGS) -o $@ $< -L$(PREFIX)/lib -l$(GMRFLibNAME) $(EXTLIBS)

example-blockupdate-with-error-control : example-blockupdate-with-error-control.o
	$(LD) $(FLAGS) -o $@ $< -L$(PREFIX)/lib -l$(GMRFLibNAME) $(EXTLIBS)

example-graph1 : example-graph1.o
	$(LD) $(FLAGS) -o $@ $< -L$(PREFIX)/lib -l$(GMRFLibNAME) $(EXTLIBS)

example-graph2 : example-graph2.o
	$(LD) $(FLAGS) -o $@ $< -L$(PREFIX)/lib -l$(GMRFLibNAME) $(EXTLIBS)

example-sample : example-sample.o
	$(LD) $(FLAGS) -o $@ $< -L$(PREFIX)/lib -l$(GMRFLibNAME) $(EXTLIBS)

example-wa : example-wa.o
	$(LD) $(FLAGS) -o $@ $< -L$(PREFIX)/lib -l$(GMRFLibNAME) $(EXTLIBS)

example-rw : example-rw.o
	$(LD) $(FLAGS) -o $@ $< -L$(PREFIX)/lib -l$(GMRFLibNAME) $(EXTLIBS)

example-geo : example-geo.o
	$(LD) $(FLAGS) -o $@ $< -L$(PREFIX)/lib -l$(GMRFLibNAME)-geo -l$(GMRFLibNAME) $(EXTLIBS)

example-qinv : example-qinv.o
	$(LD) $(FLAGS) -o $@ $< -L$(PREFIX)/lib -l$(GMRFLibNAME) $(EXTLIBS)

example-sphere : example-sphere.o
	$(LD) $(FLAGS) -o $@ $< -L$(PREFIX)/lib -l$(GMRFLibNAME) $(EXTLIBS)

example-error-handler : example-error-handler.o
	$(LD) $(FLAGS) -o $@ $< -L$(PREFIX)/lib -l$(GMRFLibNAME) $(EXTLIBS)

example-approx-1 : example-approx-1.o
	$(LD) $(FLAGS) -o $@ $< -L$(PREFIX)/lib -l$(GMRFLibNAME) $(EXTLIBS)

example-approx-2 : example-approx-2.o
	$(LD) $(FLAGS) -o $@ $< -L$(PREFIX)/lib -l$(GMRFLibNAME) $(EXTLIBS)

example-approx-3 : example-approx-3.o
	$(LD) $(FLAGS) -o $@ $< -L$(PREFIX)/lib -l$(GMRFLibNAME) $(EXTLIBS)

example-hgmrfm-1 : example-hgmrfm-1.o
	$(LD) $(FLAGS) -o $@ $< -L$(PREFIX)/lib -l$(GMRFLibNAME) $(EXTLIBS)

example-hgmrfm-2 : example-hgmrfm-2.o
	$(LD) $(FLAGS) -o $@ $< -L$(PREFIX)/lib -l$(GMRFLibNAME) $(EXTLIBS)

example-blockupdate-2 : example-blockupdate-2.o
	$(LD) $(FLAGS) -o $@ $< -L$(PREFIX)/lib -l$(GMRFLibNAME) $(EXTLIBS)

example-matern2d : example-matern2d.o
	$(LD) $(FLAGS) -o $@ $< -L$(PREFIX)/lib -l$(GMRFLibNAME) $(EXTLIBS)

../example-doxygen-%.txt : example-%.c
	expand $< > $@

doc: $(EXAMPLES_MANUAL)

clean :; rm $(PROGS) *.o sample.dat Q.dat variances.dat trace-g0.dat trace-g1.dat trace-mu.dat trace-prec.dat matern2d-sample.dat matern2d-variance.dat core.???*

.PHONY: doc clean

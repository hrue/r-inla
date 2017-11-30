SRC = bd0.c  beta.c  chebyshev.c  cospi.c  dbeta.c  dbinom.c  fmax2.c  gamma.c  lbeta.c  lgamma.c  lgammacor.c  stirlerr.c qbeta.c pbeta.c toms708.c fmin2.c pgamma.c dpois.c dnorm.c pnorm.c i1mach.c d1mach.c

from-R.c : $(SRC)
	cat $(SRC) > $@

test :; gcc  -I/usr/include/R -I. testme.c from-R.c  -lm -lgsl -lgslcblas -lm; ./a.out

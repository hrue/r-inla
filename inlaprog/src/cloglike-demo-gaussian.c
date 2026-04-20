#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <strings.h>

#include "cgeneric.h"

#define Calloc(n_, type_)  (type_ *)calloc((n_), sizeof(type_))
#define Malloc(n_, type_)  (type_ *)malloc((n_) * sizeof(type_))
#define SQR(x_) ((x_)*(x_))

double *inla_cloglike_gaussian(inla_cloglike_cmd_tp cmd, double *theta,
			       inla_cgeneric_data_tp *data, int ny, double *y, int nx, double *x, double *result)
{
#define LOG_NORMC_GAUSSIAN (-0.91893853320467274178032973640560)	/* -1/2 * log(2*pi) */
#define CDF(x_) (0.5 * (1.0 + erf(M_SQRT1_2 * (x_))))

	double *ret = NULL, prec, lprec;

	if (theta) {
		lprec = theta[0];
		prec = exp(lprec);
	} else {
		prec = lprec = NAN;
	}

	switch (cmd) {
	case INLA_CLOGLIKE_INITIAL:
	{
		// return c(M, initials)
		// where M is the number of hyperparameters

		ret = Malloc(2, double);
		ret[0] = 1;
		ret[1] = 4.0;
	}
		break;

	case INLA_CLOGLIKE_LOG_PRIOR:
	{
		// return c(LOG_PRIOR). with a Gamma(1,1) for precision, this is the log prior for the log(precision).
		ret = Malloc(1, double);
		ret[0] = -prec + lprec;
	}
		break;

	case INLA_CLOGLIKE_LOGLIKE:
	{
		// y[0] ~ N(x[i], prec)
#pragma omp simd
		for (int i = 0; i < nx; i++) {
			result[i] = LOG_NORMC_GAUSSIAN + 0.5 * (lprec - (SQR(y[0] - x[i]) * prec));
		}
	}
		break;

	case INLA_CLOGLIKE_CDF:
	{
		// Prob(y[0] < x[i]) when y[0] ~ N(x[i], prec)
		double sprec = sqrt(prec);
		for (int i = 0; i < nx; i++) {
			double z = (y[0] - x[i]) * sprec;
			result[i] = CDF(z);
		}
	}
		break;

	case INLA_CLOGLIKE_QUIT:
		break;
	}

#undef LOG_NORMC_GAUSSIAN
#undef CDF
	return (ret);
}

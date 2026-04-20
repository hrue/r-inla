#include <math.h>
#include <stdlib.h>
#include "cgeneric.h"

#define Calloc(n_, type_)  (type_ *)calloc((n_), sizeof(type_))
#define Malloc(n_, type_)  (type_ *)malloc((n_) * sizeof(type_))
#define SQR(x_) ((x_)*(x_))

double *inla_cloglike_beta(inla_cloglike_cmd_tp cmd, double *theta,
			   inla_cgeneric_data_tp *data, int ny, double *y, int nx, double *x, double *result)
{
	double *ret = NULL;

	switch (cmd) {
	case INLA_CLOGLIKE_INITIAL:
	{
		ret = Malloc(3, double);
		ret[0] = 2;
		ret[1] = 1;
		ret[2] = 1;
	}
		break;

	case INLA_CLOGLIKE_LOG_PRIOR:
	{
		double prec[2] = { 1, 1 };
		ret = Malloc(1, double);
		ret[0] = log(1.0 / sqrt(2.0 * M_PI)) + 0.5 * log(prec[0])
		    - 0.5 * prec[0] * SQR(theta[0]);
		ret[0] += log(1.0 / sqrt(2.0 * M_PI)) + 0.5 * log(prec[1])
		    - 0.5 * prec[1] * SQR(theta[1]);
	}
		break;

	case INLA_CLOGLIKE_LOGLIKE:
	{
		// covariate for phi is y[1]
		for (int i = 0; i < nx; i++) {
			double mu = exp(x[i]) / (1.0 + exp(x[i]));
			double phi = exp(theta[0] + theta[1] * y[1]);
			double a = mu * phi;
			double b = -mu * phi + phi;
			result[i] = (a - 1.0) * log(y[0])
			    + (b - 1.0) * log(1.0 - y[0])
			    + (lgamma(a + b) - (lgamma(a) + lgamma(b)));
		}
	}
		break;

	case INLA_CLOGLIKE_CDF:
	{
		for (int i = 0; i < nx; i++) {
			result[i] = NAN;
		}
	}
		break;

	case INLA_CLOGLIKE_QUIT:
		break;
	}

	return (ret);
}

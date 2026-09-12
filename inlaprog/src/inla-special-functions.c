double inla_cdf_normal(double x)
{
	/*
	 * the un-log version of inla_logcdf_normal 
	 */
	if (ABS(x) < 7.0) {
		return GMRFLib_cdfnorm(x);
	} else {
		return exp(inla_logcdf_normal(x));
	}
}

double inla_logitcdf_normal(double x)
{
	// return log(Phi(x)/(1-Phi(x)))

	if (ABS(x) < 7.0) {
		double y = inla_cdf_normal(x);
		return (log(y / (1.0 - y)));
	} else {
		// > asympt(log(Phi(x)/(1-Phi(x))), x, 16); 
		// 2
		// x 1/2 1/2 1
		// ---- + ln(x) + ln(2 Pi ) + O(----)
		// 2 2
		// 

		double val = (SQR(x) / 2.0 + log(x) + M_LN_SQRT_2PI);
		return (x > 0.0 ? val : -val);
	}
}

double inla_logcdf_normal(double x)
{
	// return the log of the cummulative distribution function for a standard normal.
	// This version is ok for all x 
	if (ABS(x) <= 7.0) {
		return (log(GMRFLib_cdfnorm(x)));
	} else {
		double t1, t4, t3, t8, t9, t13, t27, t28, t31, t47;

		if (x > 7.0) {
			t1 = 1.77245385090551602729816748334;
			t3 = M_SQRT2;
			t4 = t3 / t1;
			t8 = x * x;
			t9 = t8 * x;
			t13 = t8 * t8;
			t27 = exp(t8);
			t28 = sqrt(t27);
			t31 = 0.1e1 / M_PI;
			t47 =
			    0.1e1 / t28 * (-0.1e1 / x * t4 / 0.2e1 + 0.1e1 / t9 * t4 / 0.2e1 - 0.3e1 / 0.2e1 / t13 / x * t4 +
					   0.15e2 / 0.2e1 / t13 / t9 * t4)
			    + 0.1e1 / t27 * (-0.1e1 / t8 * t31 / 0.4e1 + 0.1e1 / t13 * t31 / 0.2e1 - 0.7e1 / 0.4e1 / t13 / t8 * t31);
			return t47;
		} else {
			// x < -7.0
			double xx = -x, cg1;
			cg1 =
			    -(pow(xx, 0.6e1) + log(0.2e1) * pow(xx, 0.4e1) + log(0.3141592653589793e1) * pow(xx, 0.4e1) +
			      0.2e1 * log(xx) * pow(xx, 0.4e1)
			      - 0.5e1 + 0.2e1 * xx * xx) * pow(xx, -0.4e1) / 0.2e1;
			return (cg1);
		}
	}
	abort();
	return 0;
}

double inla_cdf_normal_fast(double x)
{
	// a faster approximation, see misc/doc/doc/approximate-cdf-normal.pdf
	if (ABS(x) <= 7.0) {
		// see misc/doc/doc/approximate-cdf-normal.pdf
		// sqrt(M_PI / 8.0) = 0.6266570686577502....
		if (x > 0.0) {
			return (0.5 + 0.5 * sqrt(ONE_mexp(-0.6266570686577502 * SQR(x))));
		} else {
			return (1.0 - (0.5 + 0.5 * sqrt(ONE_mexp(-0.6266570686577502 * SQR(x)))));
		}
		abort();
		return (0.5 + 0.5 * sqrt(ONE_mexp(-0.6266570686577502 * SQR(x))));
	} else {
		return inla_cdf_normal(x);
	}
}

double inla_logcdf_normal_fast(double x)
{
	// a faster approximation, see misc/doc/doc/approximate-cdf-normal.pdf
	// sqrt(M_PI / 8.0) = 0.6266570686577502....
	// log(1.0/4.0) = -1.386294361119891...
	if (ABS(x) < 7.0) {
		return (log(inla_cdf_normal_fast(x)));
	} else {
		if (x > 7.0) {
			return (-0.25 * exp(-0.6266570686577502 * SQR(x)));
		} else {
			// return (log(1.0 / 4.0) - 0.6266570686577502 * SQR(x));
			return (-1.386294361119891 - 0.6266570686577502 * SQR(x));
		}
	}
}

double inla_lgamma_fast1(double x)
{
	// this is the G.Nemes (2007) approximation from https://en.wikipedia.org/wiki/Stirling's_approximation

	if (round(x) == x) {
		return gsl_sf_lnfact((int) x - 1);
	}

	double val;
	if (x < 1.0) {
		val = LGAMMAfn(x);
	} else {
		double lx = log(x);
		val = 0.5 * (LOG2PI - lx) + x * (log(x + 1.0 / (12.0 * x - 0.1 / x)) - 1.0);
	}
	return (val);
}

double inla_gamma_fast1(double x)
{
	return (exp(inla_lgamma_fast1(x)));
}

double inla_lgamma_fast2(double x)
{
	if (unlikely(x <= 0.0)) {
		return lgamma(x);
	}
#define G 5
#define N 7
	static const double p[N] = {
		1.0000000001900148240,
		76.180091729471463483,
		-86.505320329416767652,
		24.014098240830910490,
		-1.2317395724501553875,
		0.0012086509738661785061,
		-5.3952393849531283785e-6
	};

	double tmp = x + G + 0.5;
	tmp -= (x + 0.5) * log(tmp);

	double ser = p[0];
	for (int i = 1; i < N; ++i) {
		ser += p[i] / (x + i);
	}
#undef G
#undef N
	return -tmp + log(2.5066282746310005 * ser / x);
}

void inla_lgamma_fast2_m(size_t m, double *restrict x, double *restrict res)
{
	// evaluate M calls to lgamma_fast2 together, assume all x[] > 0. this is what is used for the lbeta(a,b) function, for
	// which all arguments are positive: lbeta(a,b) := lgamma(a)+lgamma(b)-lgamma(a+b)
#define G 5
#define N 7
	static const double p[N] = {
		1.0000000001900148240,
		76.180091729471463483,
		-86.505320329416767652,
		24.014098240830910490,
		-1.2317395724501553875,
		0.0012086509738661785061,
		-5.3952393849531283785e-6
	};
	double tmp[m];
	for (size_t i = 0; i < m; i++) {
		tmp[i] = x[i] + G + 0.5;
		tmp[i] -= (x[i] + 0.5) * log(tmp[i]);
	}
	double ser[m];
	GMRFLib_dfill((int) m, p[0], ser);
	for (size_t j = 1; j < N; j++) {
		for (size_t i = 0; i < m; i++) {
			ser[i] += p[j] / (x[i] + j);
		}
	}
	for (size_t i = 0; i < m; i++) {
		res[i] = -tmp[i] + log(2.5066282746310005 * ser[i] / x[i]);
	}
#undef G
#undef N
}

double inla_gamma_fast2(double x)
{
	return (exp(inla_lgamma_fast2(x)));
}


double inla_lbeta(double a, double b)
{
	double x[3] = { a, b, a + b };
	double res[3] = { 0 };
	inla_lgamma_fast2_m(3, x, res);
	return res[0] + res[1] - res[2];
}

double inla_beta(double a, double b)
{
	return exp(inla_lbeta(a, b));
}

void inla_lbeta_m(size_t m, double *restrict a, double *restrict b, double *restrict llbeta)
{
	// this is where lgamma_fast2_m() is used
	double x[3 * m];
	for (size_t i = 0, j = 0; i < 3 * m; i += 3, j++) {
		x[i] = a[j];
		x[i + 1] = b[j];
		x[i + 2] = a[j] + b[j];
	}
	double r[3 * m];
	inla_lgamma_fast2_m(3 * m, x, r);
	for (size_t i = 0, j = 0; i < 3 * m; i += 3, j++) {
		llbeta[j] = r[i] + r[i + 1] - r[i + 2];
	}
}

double inla_ipow(double x, int k)
{
	// x^k
	return gsl_sf_pow_int(x, k);
}

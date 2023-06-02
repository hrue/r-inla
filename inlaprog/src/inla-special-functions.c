
/* inla-special-functions.c
 * 
 * Copyright (C) 2007-2023 Havard Rue
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 *  * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * The author's contact information:
 *
 *        Haavard Rue
 *        CEMSE Division
 *        King Abdullah University of Science and Technology
 *        Thuwal 23955-6900, Saudi Arabia
 *        Email: haavard.rue@kaust.edu.sa
 *        Office: +966 (0)12 808 0640
 *
 */


double inla_Phi(double x)
{
	/*
	 * the un-log version of inla_log_Phi 
	 */
	if (ABS(x) < 7.0) {
		return GMRFLib_cdfnorm(x);
	} else {
		return exp(inla_log_Phi(x));
	}
}

double inla_logit_Phi(double x)
{
	// return log(Phi(x)/(1-Phi(x)))

	if (ABS(x) < 7.0) {
		double y = inla_Phi(x);
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

double inla_log_Phi(double x)
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

double inla_Phi_fast(double x)
{
	// a faster approximation, see misc/doc/doc/approximate-cdf-normal.pdf
	if (ABS(x) <= 7.0) {
		// see misc/doc/doc/approximate-cdf-normal.pdf
		// sqrt(M_PI / 8.0) = 0.6266570686577502....
		if (x > 0.0) {
			return (0.5 + 0.5 * sqrt(ONE_MINUS_EXP(-0.6266570686577502 * SQR(x))));
		} else {
			return (1.0 - (0.5 + 0.5 * sqrt(ONE_MINUS_EXP(-0.6266570686577502 * SQR(x)))));
		}
		abort();
		return (0.5 + 0.5 * sqrt(ONE_MINUS_EXP(-0.6266570686577502 * SQR(x))));
	} else {
		return inla_Phi(x);
	}
}

double inla_log_Phi_fast(double x)
{
	// a faster approximation, see misc/doc/doc/approximate-cdf-normal.pdf
	// sqrt(M_PI / 8.0) = 0.6266570686577502....
	// log(1.0/4.0) = -1.386294361119891...
	if (ABS(x) < 7.0) {
		return (log(inla_Phi_fast(x)));
	} else {
		if (x > 7.0) {
			return (-0.25 * exp(-0.6266570686577502 * SQR(x)));
		} else {
			// return (log(1.0 / 4.0) - 0.6266570686577502 * SQR(x));
			return (-1.386294361119891 - 0.6266570686577502 * SQR(x));
		}
	}
}

double inla_lgamma_fast(double x)
{
	// this is the Gergo Nemes (2007) approximation from
	// https://en.wikipedia.org/wiki/Stirling's_approximation

	if (round(x) == x) {
		return (gsl_sf_lngamma(x));
	}

	double val;
	if (x < 1.0) {
		val = gsl_sf_lngamma(x);
	} else {
		double lx = log(x), l2pi = 1.837877066409345;  // ln 2\pi
		val = 0.5 * (l2pi - lx) + x * (log(x + 1.0 / (12.0 * x - 0.1 / x)) - 1.0);
	}
	return (val);
}

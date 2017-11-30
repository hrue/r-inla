/*
 *  AUTHOR
 *	Catherine Loader, catherine@research.bell-labs.com.
 *	October 23, 2000.
 *
 *  Merge in to R:
 *	Copyright (C) 2000-2014 The R Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 *
 *
 *  DESCRIPTION
 *	Evaluates the "deviance part"
 *	bd0(x,M) :=  M * D0(x/M) = M*[ x/M * log(x/M) + 1 - (x/M) ] =
 *		  =  x * log(x/M) + M - x
 *	where M = E[X] = n*p (or = lambda), for	  x, M > 0
 *
 *	in a manner that should be stable (with small relative error)
 *	for all x and M=np. In particular for x/np close to 1, direct
 *	evaluation fails, and evaluation is based on the Taylor series
 *	of log((1+v)/(1-v)) with v = (x-M)/(x+M) = (x-np)/(x+np).
 */
#include "from-R.h"

double attribute_hidden bd0(double x, double np)
{
    double ej, s, s1, v;
    int j;

    if(!R_FINITE(x) || !R_FINITE(np) || np == 0.0) ML_ERR_return_NAN;

    if (fabs(x-np) < 0.1*(x+np)) {
	v = (x-np)/(x+np);  // might underflow to 0
	s = (x-np)*v;/* s using v -- change by MM */
	if(fabs(s) < DBL_MIN) return s;
	ej = 2*x*v;
	v = v*v;
	for (j = 1; j < 1000; j++) { /* Taylor series; 1000: no infinite loop
					as |v| < .1,  v^2000 is "zero" */
	    ej *= v;// = v^(2j+1)
	    s1 = s+ej/((j<<1)+1);
	    if (s1 == s) /* last term was effectively 0 */
		return s1 ;
	    s = s1;
	}
    }
    /* else:  | x - np |  is not too small */
    return(x*log(x/np)+np-x);
}
/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998 Ross Ihaka
 *  Copyright (C) 2000-2014 The R Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 *
 *  SYNOPSIS
 *
 *    #include <Rmath.h>
 *    double beta(double a, double b);
 *
 *  DESCRIPTION
 *
 *    This function returns the value of the beta function
 *    evaluated with arguments a and b.
 *
 *  NOTES
 *
 *    This routine is a translation into C of a Fortran subroutine
 *    by W. Fullerton of Los Alamos Scientific Laboratory.
 *    Some modifications have been made so that the routines
 *    conform to the IEEE 754 standard.
 */

#include <float.h>
#include <math.h>
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif
#include <stdlib.h>

#include "from-R.h"

double beta(double a, double b)
{
/* For IEEE double precision DBL_EPSILON = 2^-52 = 2.220446049250313e-16 :
 *   xmin, xmax : see ./gammalims.c
 *   lnsml = log(DBL_MIN) = log(2 ^ -1022) = -1022 * log(2)
 */
# define xmin  -170.5674972726612
# define xmax   171.61447887182298
# define lnsml -708.39641853226412



#ifdef IEEE_754
	/* NaNs propagated correctly */
	if(ISNAN(a) || ISNAN(b)) return a + b;
#endif


	if (a + b < xmax) {/* ~= 171.61 for IEEE */
//	return gammafn(a) * gammafn(b) / gammafn(a+b);
		/* All the terms are positive, and all can be large for large
		   or small arguments.  They are never much less than one.
		   gammafn(x) can still overflow for x ~ 1e-308,
		   but the result would too.
		*/
		return (1 / gammafn(a+b)) * gammafn(a) * gammafn(b);
	} else {
		double val = lbeta(a, b);
// underflow to 0 is not harmful per se;  exp(-999) also gives no warning
#ifndef IEEE_754
		if (val < lnsml) {
			/* a and/or b so big that beta underflows */
			ML_ERROR(ME_UNDERFLOW, "beta");
			/* return ML_UNDERFLOW; pointless giving incorrect value */
		}
#endif
		return exp(val);
	}
}
# undef xmin 
# undef xmax 
# undef lnsml 
/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998 Ross Ihaka
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 *
 *  SYNOPSIS
 *
 *    int chebyshev_init(double *dos, int nos, double eta)
 *    double chebyshev_eval(double x, double *a, int n)
 *
 *  DESCRIPTION
 *
 *    "chebyshev_init" determines the number of terms for the
 *    double precision orthogonal series "dos" needed to insure
 *    the error is no larger than "eta".  Ordinarily eta will be
 *    chosen to be one-tenth machine precision.
 *
 *    "chebyshev_eval" evaluates the n-term Chebyshev series
 *    "a" at "x".
 *
 *  NOTES
 *
 *    These routines are translations into C of Fortran routines
 *    by W. Fullerton of Los Alamos Scientific Laboratory.
 *
 *    Based on the Fortran routine dcsevl by W. Fullerton.
 *    Adapted from R. Broucke, Algorithm 446, CACM., 16, 254 (1973).
 */

#include "from-R.h"

/* NaNs propagated correctly */


int attribute_hidden chebyshev_init(double *dos, int nos, double eta)
{
    int i, ii;
    double err;

    if (nos < 1)
	return 0;

    err = 0.0;
    i = 0;			/* just to avoid compiler warnings */
    for (ii=1; ii<=nos; ii++) {
	i = nos - ii;
	err += fabs(dos[i]);
	if (err > eta) {
	    return i;
	}
    }
    return i;
}


double attribute_hidden chebyshev_eval(double x, const double *a, const int n)
{
    double b0, b1, b2, twox;
    int i;

    if (n < 1 || n > 1000) ML_ERR_return_NAN;

    if (x < -1.1 || x > 1.1) ML_ERR_return_NAN;

    twox = x * 2;
    b2 = b1 = 0;
    b0 = 0;
    for (i = 1; i <= n; i++) {
	b2 = b1;
	b1 = b0;
	b0 = twox * b1 - b2 + a[n - i];
    }
    return (b0 - b2) * 0.5;
}
/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 2013-2016 The R Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 */

#include "from-R.h"

/* HAVE_COSPI etc will not be defined in standalone-use: the
   intention is to make the versions here available in that case.

   The __cospi etc variants are from macOS (and perhaps other BSD-based systems).
*/

#ifdef HAVE_COSPI
#elif defined HAVE___COSPI
double cospi(double x) {
    return __cospi(x);
}
#else
// cos(pi * x)  -- exact when x = k/2  for all integer k
double cospi(double x) {
#ifdef IEEE_754
    /* NaNs propagated correctly */
    if (ISNAN(x)) return x;
#endif
    if(!R_FINITE(x)) ML_ERR_return_NAN;

    x = fmod(fabs(x), 2.);// cos() symmetric; cos(pi(x + 2k)) == cos(pi x) for all integer k
    if(fmod(x, 1.) == 0.5) return 0.;
    if( x == 1.)	return -1.;
    if( x == 0.)	return  1.;
    // otherwise
    return cos(M_PI * x);
}
#endif

#ifdef HAVE_SINPI
#elif defined HAVE___SINPI
double sinpi(double x) {
    return __sinpi(x);
}
#else
// sin(pi * x)  -- exact when x = k/2  for all integer k
double sinpi(double x) {
#ifdef IEEE_754
    if (ISNAN(x)) return x;
#endif
    if(!R_FINITE(x)) ML_ERR_return_NAN;

    x = fmod(x, 2.); // sin(pi(x + 2k)) == sin(pi x)  for all integer k
    // map (-2,2) --> (-1,1] :
    if(x <= -1) x += 2.; else if (x > 1.) x -= 2.;
    if(x == 0. || x == 1.) return 0.;
    if(x ==  0.5)	return  1.;
    if(x == -0.5)	return -1.;
    // otherwise
    return sin(M_PI * x);
}
#endif

// tan(pi * x)  -- exact when x = k/2  for all integer k
#if defined(HAVE_TANPI) || defined(HAVE___TANPI)
// for use in arithmetic.c, half-values documented to give NaN
double Rtanpi(double x)
#else
double tanpi(double x)
#endif
{
#ifdef IEEE_754
    if (ISNAN(x)) return x;
#endif
    if(!R_FINITE(x)) ML_ERR_return_NAN;

    x = fmod(x, 1.); // tan(pi(x + k)) == tan(pi x)  for all integer k
    // map (-1,1) --> (-1/2, 1/2] :
    if(x <= -0.5) x++; else if(x > 0.5) x--;
    return (x == 0.) ? 0. : ((x == 0.5) ? ML_NAN : tan(M_PI * x));
}

#if !defined(HAVE_TANPI) && defined(HAVE___TANPI)
double tanpi(double x) {
    return __tanpi(x);
}
#endif
/*
 *  AUTHOR
 *    Catherine Loader, catherine@research.bell-labs.com.
 *    October 23, 2000.
 *
 *  Merge in to R:
 *	Copyright (C) 2000, The R Core Team
 *  Changes to case a, b < 2, use logs to avoid underflow
 *	Copyright (C) 2006-2014 The R Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 *
 *
 *  DESCRIPTION
 *    Beta density,
 *                   (a+b-1)!     a-1       b-1
 *      p(x;a,b) = ------------ x     (1-x)
 *                 (a-1)!(b-1)!
 *
 *               = (a+b-1) dbinom(a-1; a+b-2,x)
 *
 *    The basic formula for the log density is thus
 *    (a-1) log x + (b-1) log (1-x) - lbeta(a, b)
 *    If either a or b <= 2 then 0 < lbeta(a, b) < 710 and so no
 *    term is large.  We use Loader's code only if both a and b > 2.
 */

#include "from-R.h"

double dbeta(double x, double a, double b, int give_log)
{
#ifdef IEEE_754
	/* NaNs propagated correctly */
	if (ISNAN(x) || ISNAN(a) || ISNAN(b)) return x + a + b;
#endif

	if (a < 0 || b < 0) ML_ERR_return_NAN;
	if (x < 0 || x > 1) return(R_D__0);

	// limit cases for (a,b), leading to point masses
	if(a == 0 || b == 0 || !R_FINITE(a) || !R_FINITE(b)) {
		if(a == 0 && b == 0) { // point mass 1/2 at each of {0,1} :
			if (x == 0 || x == 1) return(ML_POSINF); else return(R_D__0);
		}
		if (a == 0 || a/b == 0) { // point mass 1 at 0
			if (x == 0) return(ML_POSINF); else return(R_D__0);
		}
		if (b == 0 || b/a == 0) { // point mass 1 at 1
			if (x == 1) return(ML_POSINF); else return(R_D__0);
		}
		// else, remaining case:  a = b = Inf : point mass 1 at 1/2
		if (x == 0.5) return(ML_POSINF); else return(R_D__0);
	}

	if (x == 0) {
		if(a > 1) return(R_D__0);
		if(a < 1) return(ML_POSINF);
		/* a == 1 : */ return(R_D_val(b));
	}
	if (x == 1) {
		if(b > 1) return(R_D__0);
		if(b < 1) return(ML_POSINF);
		/* b == 1 : */ return(R_D_val(a));
	}

	double lval;
	if (a <= 2 || b <= 2)
		lval = (a-1)*log(x) + (b-1)*log1p(-x) - lbeta(a, b);
	else
		lval = log(a+b-1) + dbinom_raw(a-1, a+b-2, x, 1-x, TRUE);

	return R_D_exp(lval);
}
/*
 * AUTHOR
 *   Catherine Loader, catherine@research.bell-labs.com.
 *   October 23, 2000.
 *
 *  Merge in to R and further tweaks :
 *	Copyright (C) 2000-2015 The R Core Team
 *	Copyright (C) 2008 The R Foundation
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 *
 *
 * DESCRIPTION
 *
 *   To compute the binomial probability, call dbinom(x,n,p).
 *   This checks for argument validity, and calls dbinom_raw().
 *
 *   dbinom_raw() does the actual computation; note this is called by
 *   other functions in addition to dbinom().
 *     (1) dbinom_raw() has both p and q arguments, when one may be represented
 *         more accurately than the other (in particular, in df()).
 *     (2) dbinom_raw() does NOT check that inputs x and n are integers. This
 *         should be done in the calling function, where necessary.
 *         -- but is not the case at all when called e.g., from df() or dbeta() !
 *     (3) Also does not check for 0 <= p <= 1 and 0 <= q <= 1 or NaN's.
 *         Do this in the calling function.
 */

#include "from-R.h"

double dbinom_raw(double x, double n, double p, double q, int give_log)
{
    double lf, lc;

    if (p == 0) return((x == 0) ? R_D__1 : R_D__0);
    if (q == 0) return((x == n) ? R_D__1 : R_D__0);

    if (x == 0) {
	if(n == 0) return R_D__1;
	lc = (p < 0.1) ? -bd0(n,n*q) - n*p : n*log(q);
	return( R_D_exp(lc) );
    }
    if (x == n) {
	lc = (q < 0.1) ? -bd0(n,n*p) - n*q : n*log(p);
	return( R_D_exp(lc) );
    }
    if (x < 0 || x > n) return( R_D__0 );

    /* n*p or n*q can underflow to zero if n and p or q are small.  This
       used to occur in dbeta, and gives NaN as from R 2.3.0.  */
    lc = stirlerr(n) - stirlerr(x) - stirlerr(n-x) - bd0(x,n*p) - bd0(n-x,n*q);

    /* f = (M_2PI*x*(n-x))/n; could overflow or underflow */
    /* Upto R 2.7.1:
     * lf = log(M_2PI) + log(x) + log(n-x) - log(n);
     * -- following is much better for  x << n : */
    lf = M_LN_2PI + log(x) + log1p(- x/n);

    return R_D_exp(lc - 0.5*lf);
}

double dbinom(double x, double n, double p, int give_log)
{
#ifdef IEEE_754
    /* NaNs propagated correctly */
    if (ISNAN(x) || ISNAN(n) || ISNAN(p)) return x + n + p;
#endif

    if (p < 0 || p > 1 || R_D_negInonint(n))
	ML_ERR_return_NAN;
    R_D_nonint_check(x);
    if (x < 0 || !R_FINITE(x)) return R_D__0;

    n = R_forceint(n);
    x = R_forceint(x);

    return dbinom_raw(x, n, p, 1-p, give_log);
}
/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 */

#include "from-R.h"

double fmax2(double x, double y)
{
#ifdef IEEE_754
	if (ISNAN(x) || ISNAN(y))
		return x + y;
#endif
	return (x < y) ? y : x;
}
/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998 Ross Ihaka
 *  Copyright (C) 2000-2013 The R Core Team
 *  Copyright (C) 2002-2004 The R Foundation
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 *
 *  SYNOPSIS
 *
 *    #include <Rmath.h>
 *    double gammafn(double x);
 *
 *  DESCRIPTION
 *
 *    This function computes the value of the gamma function.
 *
 *  NOTES
 *
 *    This function is a translation into C of a Fortran subroutine
 *    by W. Fullerton of Los Alamos Scientific Laboratory.
 *    (e.g. http://www.netlib.org/slatec/fnlib/gamma.f)
 *
 *    The accuracy of this routine compares (very) favourably
 *    with those of the Sun Microsystems portable mathematical
 *    library.
 *
 *    MM specialized the case of  n!  for n < 50 - for even better precision
 */

#include "from-R.h"

double gammafn(double x)
{
    const static double gamcs[42] = {
	+.8571195590989331421920062399942e-2,
	+.4415381324841006757191315771652e-2,
	+.5685043681599363378632664588789e-1,
	-.4219835396418560501012500186624e-2,
	+.1326808181212460220584006796352e-2,
	-.1893024529798880432523947023886e-3,
	+.3606925327441245256578082217225e-4,
	-.6056761904460864218485548290365e-5,
	+.1055829546302283344731823509093e-5,
	-.1811967365542384048291855891166e-6,
	+.3117724964715322277790254593169e-7,
	-.5354219639019687140874081024347e-8,
	+.9193275519859588946887786825940e-9,
	-.1577941280288339761767423273953e-9,
	+.2707980622934954543266540433089e-10,
	-.4646818653825730144081661058933e-11,
	+.7973350192007419656460767175359e-12,
	-.1368078209830916025799499172309e-12,
	+.2347319486563800657233471771688e-13,
	-.4027432614949066932766570534699e-14,
	+.6910051747372100912138336975257e-15,
	-.1185584500221992907052387126192e-15,
	+.2034148542496373955201026051932e-16,
	-.3490054341717405849274012949108e-17,
	+.5987993856485305567135051066026e-18,
	-.1027378057872228074490069778431e-18,
	+.1762702816060529824942759660748e-19,
	-.3024320653735306260958772112042e-20,
	+.5188914660218397839717833550506e-21,
	-.8902770842456576692449251601066e-22,
	+.1527474068493342602274596891306e-22,
	-.2620731256187362900257328332799e-23,
	+.4496464047830538670331046570666e-24,
	-.7714712731336877911703901525333e-25,
	+.1323635453126044036486572714666e-25,
	-.2270999412942928816702313813333e-26,
	+.3896418998003991449320816639999e-27,
	-.6685198115125953327792127999999e-28,
	+.1146998663140024384347613866666e-28,
	-.1967938586345134677295103999999e-29,
	+.3376448816585338090334890666666e-30,
	-.5793070335782135784625493333333e-31
    };

    int i, n;
    double y;
    double sinpiy, value;

#ifdef NOMORE_FOR_THREADS
    static int ngam = 0;
    static double xmin = 0, xmax = 0., xsml = 0., dxrel = 0.;

    /* Initialize machine dependent constants, the first time gamma() is called.
	FIXME for threads ! */
    if (ngam == 0) {
	ngam = chebyshev_init(gamcs, 42, DBL_EPSILON/20);/*was .1*d1mach(3)*/
	gammalims(&xmin, &xmax);/*-> ./gammalims.c */
	xsml = exp(fmax2(log(DBL_MIN), -log(DBL_MAX)) + 0.01);
	/*   = exp(.01)*DBL_MIN = 2.247e-308 for IEEE */
	dxrel = sqrt(DBL_EPSILON);/*was sqrt(d1mach(4)) */
    }
#else
/* For IEEE double precision DBL_EPSILON = 2^-52 = 2.220446049250313e-16 :
 * (xmin, xmax) are non-trivial, see ./gammalims.c
 * xsml = exp(.01)*DBL_MIN
 * dxrel = sqrt(DBL_EPSILON) = 2 ^ -26
*/
# define ngam 22
# define xmin -170.5674972726612
# define xmax  171.61447887182298
# define xsml 2.2474362225598545e-308
# define dxrel 1.490116119384765696e-8
#endif

    if(ISNAN(x)) return x;

    /* If the argument is exactly zero or a negative integer
     * then return NaN. */
    if (x == 0 || (x < 0 && x == round(x))) {
	ML_ERROR(ME_DOMAIN, "gammafn");
	return ML_NAN;
    }

    y = fabs(x);

    if (y <= 10) {

	/* Compute gamma(x) for -10 <= x <= 10
	 * Reduce the interval and find gamma(1 + y) for 0 <= y < 1
	 * first of all. */

	n = (int) x;
	if(x < 0) --n;
	y = x - n;/* n = floor(x)  ==>	y in [ 0, 1 ) */
	--n;
	value = chebyshev_eval(y * 2 - 1, gamcs, ngam) + .9375;
	if (n == 0)
	    return value;/* x = 1.dddd = 1+y */

	if (n < 0) {
	    /* compute gamma(x) for -10 <= x < 1 */

	    /* exact 0 or "-n" checked already above */

	    /* The answer is less than half precision */
	    /* because x too near a negative integer. */
	    if (x < -0.5 && fabs(x - (int)(x - 0.5) / x) < dxrel) {
		ML_ERROR(ME_PRECISION, "gammafn");
	    }

	    /* The argument is so close to 0 that the result would overflow. */
	    if (y < xsml) {
		ML_ERROR(ME_RANGE, "gammafn");
		if(x > 0) return ML_POSINF;
		else return ML_NEGINF;
	    }

	    n = -n;

	    for (i = 0; i < n; i++) {
		value /= (x + i);
	    }
	    return value;
	}
	else {
	    /* gamma(x) for 2 <= x <= 10 */

	    for (i = 1; i <= n; i++) {
		value *= (y + i);
	    }
	    return value;
	}
    }
    else {
	/* gamma(x) for	 y = |x| > 10. */

	if (x > xmax) {			/* Overflow */
	    ML_ERROR(ME_RANGE, "gammafn");
	    return ML_POSINF;
	}

	if (x < xmin) {			/* Underflow */
	    ML_ERROR(ME_UNDERFLOW, "gammafn");
	    return 0.;
	}

	if(y <= 50 && y == (int)y) { /* compute (n - 1)! */
	    value = 1.;
	    for (i = 2; i < y; i++) value *= i;
	}
	else { /* normal case */
	    value = exp((y - 0.5) * log(y) - y + M_LN_SQRT_2PI +
			((2*y == (int)2*y)? stirlerr(y) : lgammacor(y)));
	}
	if (x > 0)
	    return value;

	if (fabs((x - (int)(x - 0.5))/x) < dxrel){

	    /* The answer is less than half precision because */
	    /* the argument is too near a negative integer. */

	    ML_ERROR(ME_PRECISION, "gammafn");
	}

	sinpiy = sinpi(y);
	if (sinpiy == 0) {		/* Negative integer arg - overflow */
	    ML_ERROR(ME_RANGE, "gammafn");
	    return ML_POSINF;
	}

	return -M_PI / (y * sinpiy * value);
    }
}
# undef xmin 
# undef xmax 
# undef lnsml 
# undef dxrel
/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998 Ross Ihaka
 *  Copyright (C) 2000-12 The R Core Team
 *  Copyright (C) 2003 The R Foundation
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 *
 *  SYNOPSIS
 *
 *    #include <Rmath.h>
 *    double lbeta(double a, double b);
 *
 *  DESCRIPTION
 *
 *    This function returns the value of the log beta function.
 *
 *  NOTES
 *
 *    This routine is a translation into C of a Fortran subroutine
 *    by W. Fullerton of Los Alamos Scientific Laboratory.
 */

#include "from-R.h"

double lbeta(double a, double b)
{
    double corr, p, q;

#ifdef IEEE_754
    if(ISNAN(a) || ISNAN(b))
	return a + b;
#endif
    p = q = a;
    if(b < p) p = b;/* := min(a,b) */
    if(b > q) q = b;/* := max(a,b) */

    /* both arguments must be >= 0 */
    if (p < 0)
	ML_ERR_return_NAN
    else if (p == 0) {
	return ML_POSINF;
    }
    else if (!R_FINITE(q)) { /* q == +Inf */
	return ML_NEGINF;
    }

    if (p >= 10) {
	/* p and q are big. */
	corr = lgammacor(p) + lgammacor(q) - lgammacor(p + q);
	return log(q) * -0.5 + M_LN_SQRT_2PI + corr
		+ (p - 0.5) * log(p / (p + q)) + q * log1p(-p / (p + q));
    }
    else if (q >= 10) {
	/* p is small, but q is big. */
	corr = lgammacor(q) - lgammacor(p + q);
	return lgammafn(p) + corr + p - p * log(p + q)
		+ (q - 0.5) * log1p(-p / (p + q));
    }
    else {
	/* p and q are small: p <= q < 10. */
	/* R change for very small args */
	if (p < 1e-306) return lgamma(p) + (lgamma(q) - lgamma(p+q));
	else return log(gammafn(p) * (gammafn(q) / gammafn(p + q)));
    }
}
/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998 Ross Ihaka
 *  Copyright (C) 2000-2012 The R Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 *
 *  SYNOPSIS
 *
 *    #include <Rmath.h>
 *    double lgammafn_sign(double x, int *sgn);
 *    double lgammafn(double x);
 *
 *  DESCRIPTION
 *
 *    The function lgammafn computes log|gamma(x)|.  The function
 *    lgammafn_sign in addition assigns the sign of the gamma function
 *    to the address in the second argument if this is not NULL.
 *
 *  NOTES
 *
 *    This routine is a translation into C of a Fortran subroutine
 *    by W. Fullerton of Los Alamos Scientific Laboratory.
 *
 *    The accuracy of this routine compares (very) favourably
 *    with those of the Sun Microsystems portable mathematical
 *    library.
 */

#include "from-R.h"

double lgammafn_sign(double x, int *sgn)
{
    double ans, y, sinpiy;

#ifdef NOMORE_FOR_THREADS
    static double xmax = 0.;
    static double dxrel = 0.;

    if (xmax == 0) {/* initialize machine dependent constants _ONCE_ */
	xmax = d1mach(2)/log(d1mach(2));/* = 2.533 e305	 for IEEE double */
	dxrel = sqrt (d1mach(4));/* sqrt(Eps) ~ 1.49 e-8  for IEEE double */
    }
#else
/* For IEEE double precision DBL_EPSILON = 2^-52 = 2.220446049250313e-16 :
   xmax  = DBL_MAX / log(DBL_MAX) = 2^1024 / (1024 * log(2)) = 2^1014 / log(2)
   dxrel = sqrt(DBL_EPSILON) = 2^-26 = 5^26 * 1e-26 (is *exact* below !)
 */
#define xmax  2.5327372760800758e+305
#define dxrel 1.490116119384765625e-8
#endif

    if (sgn != NULL) *sgn = 1;

#ifdef IEEE_754
    if(ISNAN(x)) return x;
#endif

    if (sgn != NULL && x < 0 && fmod(floor(-x), 2.) == 0)
	*sgn = -1;

    if (x <= 0 && x == trunc(x)) { /* Negative integer argument */
	ML_ERROR(ME_RANGE, "lgamma");
	return ML_POSINF;/* +Inf, since lgamma(x) = log|gamma(x)| */
    }

    y = fabs(x);

    if (y < 1e-306) return -log(y); // denormalized range, R change
    if (y <= 10) return log(fabs(gammafn(x)));
    /*
      ELSE  y = |x| > 10 ---------------------- */

    if (y > xmax) {
	ML_ERROR(ME_RANGE, "lgamma");
	return ML_POSINF;
    }

    if (x > 0) { /* i.e. y = x > 10 */
#ifdef IEEE_754
	if(x > 1e17)
	    return(x*(log(x) - 1.));
	else if(x > 4934720.)
	    return(M_LN_SQRT_2PI + (x - 0.5) * log(x) - x);
	else
#endif
	    return M_LN_SQRT_2PI + (x - 0.5) * log(x) - x + lgammacor(x);
    }
    /* else: x < -10; y = -x */
    sinpiy = fabs(sinpi(y));

    if (sinpiy == 0) { /* Negative integer argument ===
			  Now UNNECESSARY: caught above */
	MATHLIB_WARNING(" ** should NEVER happen! *** [lgamma.c: Neg.int, y=%g]\n",y);
	ML_ERR_return_NAN;
    }

    ans = M_LN_SQRT_PId2 + (x - 0.5) * log(y) - x - log(sinpiy) - lgammacor(y);

    if(fabs((x - trunc(x - 0.5)) * ans / x) < dxrel) {

	/* The answer is less than half precision because
	 * the argument is too near a negative integer. */

	ML_ERROR(ME_PRECISION, "lgamma");
    }

    return ans;
}

double lgammafn(double x)
{
    return lgammafn_sign(x, NULL);
}
# undef xmin 
# undef xmax 
# undef lnsml 
# undef dxrel
/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998 Ross Ihaka
 *  Copyright (C) 2000-2001 The R Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 *
 *  SYNOPSIS
 *
 *    #include <Rmath.h>
 *    double lgammacor(double x);
 *
 *  DESCRIPTION
 *
 *    Compute the log gamma correction factor for x >= 10 so that
 *
 *    log(gamma(x)) = .5*log(2*pi) + (x-.5)*log(x) -x + lgammacor(x)
 *
 *    [ lgammacor(x) is called	Del(x)	in other contexts (e.g. dcdflib)]
 *
 *  NOTES
 *
 *    This routine is a translation into C of a Fortran subroutine
 *    written by W. Fullerton of Los Alamos Scientific Laboratory.
 *
 *  SEE ALSO
 *
 *    Loader(1999)'s stirlerr() {in ./stirlerr.c} is *very* similar in spirit,
 *    is faster and cleaner, but is only defined "fast" for half integers.
 */

#include "from-R.h"

double attribute_hidden lgammacor(double x)
{
    const static double algmcs[15] = {
	+.1666389480451863247205729650822e+0,
	-.1384948176067563840732986059135e-4,
	+.9810825646924729426157171547487e-8,
	-.1809129475572494194263306266719e-10,
	+.6221098041892605227126015543416e-13,
	-.3399615005417721944303330599666e-15,
	+.2683181998482698748957538846666e-17,
	-.2868042435334643284144622399999e-19,
	+.3962837061046434803679306666666e-21,
	-.6831888753985766870111999999999e-23,
	+.1429227355942498147573333333333e-24,
	-.3547598158101070547199999999999e-26,
	+.1025680058010470912000000000000e-27,
	-.3401102254316748799999999999999e-29,
	+.1276642195630062933333333333333e-30
    };

    double tmp;

/* For IEEE double precision DBL_EPSILON = 2^-52 = 2.220446049250313e-16 :
 *   xbig = 2 ^ 26.5
 *   xmax = DBL_MAX / 48 =  2^1020 / 3 */
#define nalgm 5
#define xbig  94906265.62425156
#define xmax  3.745194030963158e306

    if (x < 10)
	ML_ERR_return_NAN
    else if (x >= xmax) {
	ML_ERROR(ME_UNDERFLOW, "lgammacor");
	/* allow to underflow below */
    }
    else if (x < xbig) {
	tmp = 10 / x;
	return chebyshev_eval(tmp * tmp * 2 - 1, algmcs, nalgm) / x;
    }
    return 1 / (x * 12);
}
# undef xmin 
# undef xmax 
# undef lnsml 
/*
 *  AUTHOR
 *    Catherine Loader, catherine@research.bell-labs.com.
 *    October 23, 2000.
 *
 *  Merge in to R:
 *	Copyright (C) 2000, The R Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 *
 *
 *  DESCRIPTION
 *
 *    Computes the log of the error term in Stirling's formula.
 *      For n > 15, uses the series 1/12n - 1/360n^3 + ...
 *      For n <=15, integers or half-integers, uses stored values.
 *      For other n < 15, uses lgamma directly (don't use this to
 *        write lgamma!)
 *
 * Merge in to R:
 * Copyright (C) 2000, The R Core Team
 * R has lgammafn, and lgamma is not part of ISO C
 */

#include "from-R.h"

/* stirlerr(n) = log(n!) - log( sqrt(2*pi*n)*(n/e)^n )
 *             = log Gamma(n+1) - 1/2 * [log(2*pi) + log(n)] - n*[log(n) - 1]
 *             = log Gamma(n+1) - (n + 1/2) * log(n) + n - log(2*pi)/2
 *
 * see also lgammacor() in ./lgammacor.c  which computes almost the same!
 */

double attribute_hidden stirlerr(double n)
{

#define S0 0.083333333333333333333       /* 1/12 */
#define S1 0.00277777777777777777778     /* 1/360 */
#define S2 0.00079365079365079365079365  /* 1/1260 */
#define S3 0.000595238095238095238095238 /* 1/1680 */
#define S4 0.0008417508417508417508417508/* 1/1188 */

/*
  error for 0, 0.5, 1.0, 1.5, ..., 14.5, 15.0.
*/
    const static double sferr_halves[31] = {
	0.0, /* n=0 - wrong, place holder only */
	0.1534264097200273452913848,  /* 0.5 */
	0.0810614667953272582196702,  /* 1.0 */
	0.0548141210519176538961390,  /* 1.5 */
	0.0413406959554092940938221,  /* 2.0 */
	0.03316287351993628748511048, /* 2.5 */
	0.02767792568499833914878929, /* 3.0 */
	0.02374616365629749597132920, /* 3.5 */
	0.02079067210376509311152277, /* 4.0 */
	0.01848845053267318523077934, /* 4.5 */
	0.01664469118982119216319487, /* 5.0 */
	0.01513497322191737887351255, /* 5.5 */
	0.01387612882307074799874573, /* 6.0 */
	0.01281046524292022692424986, /* 6.5 */
	0.01189670994589177009505572, /* 7.0 */
	0.01110455975820691732662991, /* 7.5 */
	0.010411265261972096497478567, /* 8.0 */
	0.009799416126158803298389475, /* 8.5 */
	0.009255462182712732917728637, /* 9.0 */
	0.008768700134139385462952823, /* 9.5 */
	0.008330563433362871256469318, /* 10.0 */
	0.007934114564314020547248100, /* 10.5 */
	0.007573675487951840794972024, /* 11.0 */
	0.007244554301320383179543912, /* 11.5 */
	0.006942840107209529865664152, /* 12.0 */
	0.006665247032707682442354394, /* 12.5 */
	0.006408994188004207068439631, /* 13.0 */
	0.006171712263039457647532867, /* 13.5 */
	0.005951370112758847735624416, /* 14.0 */
	0.005746216513010115682023589, /* 14.5 */
	0.005554733551962801371038690  /* 15.0 */
    };
    double nn;

    if (n <= 15.0) {
	nn = n + n;
	if (nn == (int)nn) return(sferr_halves[(int)nn]);
	return(lgammafn(n + 1.) - (n + 0.5)*log(n) + n - M_LN_SQRT_2PI);
    }

    nn = n*n;
    if (n>500) return((S0-S1/nn)/n);
    if (n> 80) return((S0-(S1-S2/nn)/nn)/n);
    if (n> 35) return((S0-(S1-(S2-S3/nn)/nn)/nn)/n);
    /* 15 < n <= 35 : */
    return((S0-(S1-(S2-(S3-S4/nn)/nn)/nn)/nn)/n);
}
/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1998--2017  The R Core Team
 *  Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 *  based on code (C) 1979 and later Royal Statistical Society
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 *

 * Reference:
 * Cran, G. W., K. J. Martin and G. E. Thomas (1977).
 *	Remark AS R19 and Algorithm AS 109,
 *	Applied Statistics, 26(1), 111-114.
 * Remark AS R83 (v.39, 309-310) and the correction (v.40(1) p.236)
 *	have been incorporated in this version.
 */


#include <R_ext/Arith.h>
#include "from-R.h"


#ifdef DEBUG_qbeta
# define R_ifDEBUG_printf(...) REprintf(__VA_ARGS__)
#else
# define R_ifDEBUG_printf(...)
#endif

#define USE_LOG_X_CUTOFF -5.
//                       --- based on some testing; had = -10

#define n_NEWTON_FREE 4
//                   --- based on some testing; had = 10

#define MLOGICAL_NA -1
// an "NA_LOGICAL" substitute for Mathlib {only used here, for now}

//attribute_hidden
static void
qbeta_raw(double alpha, double p, double q, int lower_tail, int log_p,
	  int swap_01, double log_q_cut, int n_N, double* qb);

double qbeta(double alpha, double p, double q, int lower_tail, int log_p)
{

    /* test for admissibility of parameters */
#ifdef IEEE_754
    if (ISNAN(p) || ISNAN(q) || ISNAN(alpha))
	return p + q + alpha;
#endif
    if(p < 0. || q < 0.) ML_ERR_return_NAN;
    // allowing p==0 and q==0  <==> treat as one- or two-point mass

    double qbet[2];// = { qbeta(), 1 - qbeta() }
    qbeta_raw(alpha, p, q, lower_tail, log_p,
	      MLOGICAL_NA, USE_LOG_X_CUTOFF, n_NEWTON_FREE, qbet);
    return qbet[0];
}

static const double
#ifdef IEEE_754
// CARE: assumes subnormal numbers, i.e., no underflow at DBL_MIN:
    DBL_very_MIN  = DBL_MIN / 4.,
    DBL_log_v_MIN = M_LN2*(DBL_MIN_EXP - 2),
// Too extreme: inaccuracy in pbeta(); e.g for  qbeta(0.95, 1e-9, 20):
// -> in pbeta() --> bgrat(..... b*z == 0 underflow, hence inaccurate pbeta()
    /* DBL_very_MIN  = 0x0.0000001p-1022, // = 2^-1050 = 2^(-1022 - 28) */
    /* DBL_log_v_MIN = -1050. * M_LN2, // = log(DBL_very_MIN) */
// the most extreme -- not ok, as pbeta() then behaves strangely,
// e.g., for  qbeta(0.95, 1e-8, 20):
    /* DBL_very_MIN  = 0x0.0000000000001p-1022, // = 2^-1074 = 2^(-1022 -52) */
    /* DBL_log_v_MIN = -1074. * M_LN2, // = log(DBL_very_MIN) */

    DBL_1__eps    = 0x1.fffffffffffffp-1; // = 1 - 2^-53
#else // untested :
    DBL_1__eps    = 1 - DBL_EPSILON;     // or rather (1 - DBL_EPSILON/2) (??)
#endif

/* set the exponent of acu to -2r-2 for r digits of accuracy */
/*---- NEW ---- -- still fails for p = 1e11, q=.5*/

#define fpu 3e-308
/* acu_min:  Minimal value for accuracy 'acu' which will depend on (a,p);
	     acu_min >= fpu ! */
#define acu_min 1e-300
#define p_lo fpu
#define p_hi 1-2.22e-16

#define const1 2.30753
#define const2 0.27061
#define const3 0.99229
#define const4 0.04481

// Returns both qbeta() and its "mirror" 1-qbeta(). Useful notably when qbeta() ~= 1
attribute_hidden void
qbeta_raw(double alpha, double p, double q, int lower_tail, int log_p,
	  int swap_01, // {TRUE, NA, FALSE}: if NA, algorithm decides swap_tail
	  double log_q_cut, /* if == Inf: return log(qbeta(..));
			       otherwise, if finite: the bound for
			       switching to log(x)-scale; see use_log_x */
	  int n_N,  // number of "unconstrained" Newton steps before switching to constrained
	  double *qb) // = qb[0:1] = { qbeta(), 1 - qbeta() }
{
    Rboolean
	swap_choose = (swap_01 == MLOGICAL_NA),
	swap_tail,
	log_, give_log_q = (log_q_cut == ML_POSINF),
	use_log_x = give_log_q, // or u < log_q_cut  below
	warned = FALSE, add_N_step = TRUE;
    int i_pb, i_inn;
    double a, la, logbeta, g, h, pp, p_, qq, r, s, t, w, y = -1.;
    volatile double u, xinbta;

    // Assuming p >= 0, q >= 0  here ...

    // Deal with boundary cases here:
    if(alpha == R_DT_0) {
#define return_q_0						\
	if(give_log_q) { qb[0] = ML_NEGINF; qb[1] = 0; }	\
	else {           qb[0] = 0;         qb[1] = 1; }	\
	return

	return_q_0;
    }
    if(alpha == R_DT_1) {
#define return_q_1						\
	if(give_log_q) { qb[0] = 0; qb[1] = ML_NEGINF; }	\
	else {           qb[0] = 1; qb[1] = 0;         }	\
	return

	return_q_1;
    }

    // check alpha {*before* transformation which may all accuracy}:
    if((log_p && alpha > 0) ||
       (!log_p && (alpha < 0 || alpha > 1))) { // alpha is outside
	R_ifDEBUG_printf("qbeta(alpha=%g, %g, %g, .., log_p=%d): %s%s\n",
			 alpha, p,q, log_p, "alpha not in ",
			 log_p ? "[-Inf, 0]" : "[0,1]");
	// ML_ERR_return_NAN :
	ML_ERROR(ME_DOMAIN, "");
	qb[0] = qb[1] = ML_NAN; return;
    }

    //  p==0, q==0, p = Inf, q = Inf  <==> treat as one- or two-point mass
    if(p == 0 || q == 0 || !R_FINITE(p) || !R_FINITE(q)) {
	// We know 0 < T(alpha) < 1 : pbeta() is constant and trivial in {0, 1/2, 1}
	R_ifDEBUG_printf(
	    "qbeta(%g, %g, %g, lower_t=%d, log_p=%d): (p,q)-boundary: trivial\n",
	    alpha, p,q, lower_tail, log_p);
	if(p == 0 && q == 0) { // point mass 1/2 at each of {0,1} :
	    if(alpha < R_D_half) { return_q_0; }
	    if(alpha > R_D_half) { return_q_1; }
	    // else:  alpha == "1/2"
#define return_q_half					\
	    if(give_log_q) qb[0] = qb[1] = -M_LN2;	\
	    else	   qb[0] = qb[1] = 0.5;		\
	    return

	    return_q_half;
	} else if (p == 0 || p/q == 0) { // point mass 1 at 0 - "flipped around"
	    return_q_0;
	} else if (q == 0 || q/p == 0) { // point mass 1 at 0 - "flipped around"
	    return_q_1;
	}
	// else:  p = q = Inf : point mass 1 at 1/2
	return_q_half;
    }

    /* initialize */
    p_ = R_DT_qIv(alpha);/* lower_tail prob (in any case) */
    // Conceptually,  0 < p_ < 1  (but can be 0 or 1 because of cancellation!)
    logbeta = lbeta(p, q);

    swap_tail = (swap_choose) ? (p_ > 0.5) : swap_01;
    // change tail; default (swap_01 = NA): afterwards 0 < a <= 1/2
    if(swap_tail) { /* change tail, swap  p <-> q :*/
	a = R_DT_CIv(alpha); // = 1 - p_ < 1/2
	/* la := log(a), but without numerical cancellation: */
	la = R_DT_Clog(alpha);
	pp = q; qq = p;
    }
    else {
	a = p_;
	la = R_DT_log(alpha);
	pp = p; qq = q;
    }

    /* calculate the initial approximation */

    /* Desired accuracy for Newton iterations (below) should depend on  (a,p)
     * This is from Remark .. on AS 109, adapted.
     * However, it's not clear if this is "optimal" for IEEE double prec.

     * acu = fmax2(acu_min, pow(10., -25. - 5./(pp * pp) - 1./(a * a)));

     * NEW: 'acu' accuracy NOT for squared adjustment, but simple;
     * ---- i.e.,  "new acu" = sqrt(old acu)
     */
    double acu = fmax2(acu_min, pow(10., -13. - 2.5/(pp * pp) - 0.5/(a * a)));
    // try to catch  "extreme left tail" early
    double tx, u0 = (la + log(pp) + logbeta) / pp; // = log(x_0)
    static const double
	log_eps_c = M_LN2 * (1. - DBL_MANT_DIG);// = log(DBL_EPSILON) = -36.04..
    r = pp*(1.-qq)/(pp+1.);

    t = 0.2;
    // FIXME: Factor 0.2 is a bit arbitrary;  '1' is clearly much too much.

    R_ifDEBUG_printf(
	"qbeta(%g, %g, %g, lower_t=%d, log_p=%d):%s\n"
	"  swap_tail=%d, la=%#8g, u0=%#8g (bnd: %g (%g)) ",
	alpha, p,q, lower_tail, log_p,
	(log_p && (p_ == 0. || p_ == 1.)) ? (p_==0.?" p_=0":" p_=1") : "",
	swap_tail, la, u0,
	(t*log_eps_c - log(fabs(pp*(1.-qq)*(2.-qq)/(2.*(pp+2.)))))/2.,
	 t*log_eps_c - log(fabs(r))
	);

    if(M_LN2 * DBL_MIN_EXP < u0 && // cannot allow exp(u0) = 0 ==> exp(u1) = exp(u0) = 0
       u0 < -0.01 && // (must: u0 < 0, but too close to 0 <==> x = exp(u0) = 0.99..)
       // qq <= 2 && // <--- "arbitrary"
       // u0 <  t*log_eps_c - log(fabs(r)) &&
       u0 < (t*log_eps_c - log(fabs(pp*(1.-qq)*(2.-qq)/(2.*(pp+2.)))))/2.)
    {
// TODO: maybe jump here from below, when initial u "fails" ?
// L_tail_u:
	// MM's one-step correction (cheaper than 1 Newton!)
	r = r*exp(u0);// = r*x0
	if(r > -1.) {
	    u = u0 - log1p(r)/pp;
	    R_ifDEBUG_printf("u1-u0=%9.3g --> choosing u = u1\n", u-u0);
	} else {
	    u = u0;
	    R_ifDEBUG_printf("cannot cheaply improve u0\n");
	}
	tx = xinbta = exp(u);
	use_log_x = TRUE; // or (u < log_q_cut)  ??
	goto L_Newton;
    }

    // y := y_\alpha in AS 64 := Hastings(1955) approximation of qnorm(1 - a) :
    r = sqrt(-2 * la);
    y = r - (const1 + const2 * r) / (1. + (const3 + const4 * r) * r);

    if (pp > 1 && qq > 1) { // use  Carter(1947), see AS 109, remark '5.'
	r = (y * y - 3.) / 6.;
	s = 1. / (pp + pp - 1.);
	t = 1. / (qq + qq - 1.);
	h = 2. / (s + t);
	w = y * sqrt(h + r) / h - (t - s) * (r + 5. / 6. - 2. / (3. * h));
	R_ifDEBUG_printf("p,q > 1 => w=%g", w);
	if(w > 300) { // exp(w+w) is huge or overflows
	    t = w+w + log(qq) - log(pp); // = argument of log1pexp(.)
	    u = // log(xinbta) = - log1p(qq/pp * exp(w+w)) = -log(1 + exp(t))
		(t <= 18) ? -log1p(exp(t)) : -t - exp(-t);
	    xinbta = exp(u);
	} else {
	    xinbta = pp / (pp + qq * exp(w + w));
	    u = // log(xinbta)
		- log1p(qq/pp * exp(w+w));
	}
    } else { // use the original AS 64 proposal, Scheffé-Tukey (1944) and Wilson-Hilferty
	r = qq + qq;
	/* A slightly more stable version of  t := \chi^2_{alpha} of AS 64
	 * t = 1. / (9. * qq); t = r * R_pow_di(1. - t + y * sqrt(t), 3);  */
	t = 1. / (3. * sqrt(qq));
	t = r * R_pow_di(1. + t*(-t + y), 3);// = \chi^2_{alpha} of AS 64
	s = 4. * pp + r - 2.;// 4p + 2q - 2 = numerator of new t = (...) / chi^2
	R_ifDEBUG_printf("min(p,q) <= 1: t=%g", t);
	if (t == 0 || (t < 0. && s >= t)) { // cannot use chisq approx
	    // x0 = 1 - { (1-a)*q*B(p,q) } ^{1/q}    {AS 65}
	    // xinbta = 1. - exp((log(1-a)+ log(qq) + logbeta) / qq);
	    double l1ma;/* := log(1-a), directly from alpha (as 'la' above):
			 * FIXME: not worth it? log1p(-a) always the same ?? */
	    if(swap_tail)
		l1ma = R_DT_log(alpha);
	    else
		l1ma = R_DT_Clog(alpha);
	    R_ifDEBUG_printf(" t <= 0 : log1p(-a)=%.15g, better l1ma=%.15g\n", log1p(-a), l1ma);
	    double xx = (l1ma + log(qq) + logbeta) / qq;
	    if(xx <= 0.) {
		xinbta = -expm1(xx);
		u = R_Log1_Exp (xx);// =  log(xinbta) = log(1 - exp(...A...))
	    } else { // xx > 0 ==> 1 - e^xx < 0 .. is nonsense
		R_ifDEBUG_printf(" xx=%g > 0: xinbta:= 1-e^xx < 0\n", xx);
		xinbta = 0; u = ML_NEGINF; /// FIXME can do better?
	    }
	} else {
	    t = s / t;
	    R_ifDEBUG_printf(" t > 0 or s < t < 0:  new t = %g ( > 1 ?)\n", t);
	    if (t <= 1.) { // cannot use chisq, either
		u = (la + log(pp) + logbeta) / pp;
		xinbta = exp(u);
	    } else { // (1+x0)/(1-x0) = t,  solved for x0 :
		xinbta = 1. - 2. / (t + 1.);
		u = log1p(-2. / (t + 1.));
	    }
	}
    }

    // Problem: If initial u is completely wrong, we make a wrong decision here
    if(swap_choose &&
       (( swap_tail && u >= -exp(  log_q_cut)) || // ==> "swap back"
	(!swap_tail && u >= -exp(4*log_q_cut) && pp / qq < 1000.) // ==> "swap now"
	   )) {
	// "revert swap" -- and use_log_x
	swap_tail = !swap_tail;
	R_ifDEBUG_printf(" u = %g (e^u = xinbta = %.16g) ==> ", u, xinbta);
	if(swap_tail) { // "swap now" (much less easily)
	    a = R_DT_CIv(alpha); // needed ?
	    la = R_DT_Clog(alpha);
	    pp = q; qq = p;
	}
	else { // swap back :
	    a = p_;
	    la = R_DT_log(alpha);
	    pp = p; qq = q;
	}
	R_ifDEBUG_printf("\"%s\"; la = %g\n",
			 (swap_tail ? "swap now" : "swap back"), la);
	// we could redo computations above, but this should be stable
	u = R_Log1_Exp(u);
	xinbta = exp(u);

/* Careful: "swap now"  should not fail if
   1) the above initial xinbta is "completely wrong"
   2) The correction step can go outside (u_n > 0 ==>  e^u > 1 is illegal)
   e.g., for  qbeta(0.2066, 0.143891, 0.05)
*/
    } else R_ifDEBUG_printf("\n");

    if(!use_log_x)
	use_log_x = (u < log_q_cut);// <==> xinbta = e^u < exp(log_q_cut)
    Rboolean
	bad_u = !R_FINITE(u),
	bad_init = bad_u || xinbta > p_hi;

    R_ifDEBUG_printf(" -> u = %g, e^u = xinbta = %.16g, (Newton acu=%g%s%s%s)\n",
		     u, xinbta, acu, (bad_u ? ", ** bad u **" : ""),
		     ((bad_init && !bad_u) ? ", ** bad_init **" : ""),
		     (use_log_x ? ", on u = LOG(x) SCALE" : ""));

    double u_n = 1.; // -Wall
    tx = xinbta; // keeping "original initial x" (for now)

    if(bad_u || u < log_q_cut) {
	/* e.g.
	   qbeta(0.21, .001, 0.05)
	   try "left border" quickly, i.e.,
	   try at smallest positive number: */
	w = pbeta_raw(DBL_very_MIN, pp, qq, TRUE, log_p);
	if(w > (log_p ? la : a)) {
	    R_ifDEBUG_printf(
		" quantile is left of %g; \"convergence\"\n", DBL_very_MIN);
	    if(log_p || fabs(w - a) < fabs(0 - a)) { // DBL_very_MIN is better than 0
		tx   = DBL_very_MIN;
		u_n  = DBL_log_v_MIN;// = log(DBL_very_MIN)
	    } else {
		tx   = 0.;
		u_n  = ML_NEGINF;
	    }
	    use_log_x = log_p; add_N_step = FALSE; goto L_return;
	}
	else {
	    R_ifDEBUG_printf(" pbeta(%g, *) = %g <= %g (= %s) --> continuing\n",
			     DBL_log_v_MIN, w, (log_p ? la : a), (log_p ? "la" : "a"));
	    if(u  < DBL_log_v_MIN) {
		u = DBL_log_v_MIN;// = log(DBL_very_MIN)
		xinbta = DBL_very_MIN;
	    }
	}
    }


    /* Sometimes the approximation is negative (and == 0 is also not "ok") */
    if (bad_init && !(use_log_x && tx > 0)) {
	if(u == ML_NEGINF) {
	    R_ifDEBUG_printf("  u = -Inf;");
	    u = M_LN2 * DBL_MIN_EXP;
	    xinbta = DBL_MIN;
	} else {
	    R_ifDEBUG_printf(" bad_init: u=%g, xinbta=%g;", u,xinbta);
	    xinbta = (xinbta > 1.1) // i.e. "way off"
		? 0.5 // otherwise, keep the respective boundary:
		: ((xinbta < p_lo) ? exp(u) : p_hi);
	    if(bad_u)
		u = log(xinbta);
	    // otherwise: not changing "potentially better" u than the above
	}
	R_ifDEBUG_printf(" -> (partly)new u=%g, xinbta=%g\n", u,xinbta);
    }

L_Newton:
    /* --------------------------------------------------------------------

     * Solve for x by a modified Newton-Raphson method, using pbeta_raw()
     */
    r = 1 - pp;
    t = 1 - qq;
    double wprev = 0., prev = 1., adj = 1.; // -Wall

    if(use_log_x) { // find  log(xinbta) -- work in  u := log(x) scale
	// if(bad_init && tx > 0) xinbta = tx;// may have been better

	for (i_pb=0; i_pb < 1000; i_pb++) {
	    // using log_p == TRUE  unconditionally here
	    /* FIXME: if exp(u) = xinbta underflows to 0,
	     *  want different formula pbeta_log(u, ..) */
	    y = pbeta_raw(xinbta, pp, qq, /*lower_tail = */ TRUE, TRUE);

	    /* w := Newton step size for   L(u) = log F(e^u)  =!= 0;   u := log(x)
	     *   =  (L(.) - la) / L'(.);  L'(u)= (F'(e^u) * e^u ) / F(e^u)
	     *   =  (L(.) - la)*F(.) / {F'(e^u) * e^u } =
	     *   =  (L(.) - la) * e^L(.) * e^{-log F'(e^u) - u}
	     *   =  ( y   - la) * e^{ y - u -log F'(e^u)}
	     and  -log F'(x)= -log f(x) = - -logbeta + (1-p) log(x) + (1-q) log(1-x)
	                                = logbeta + (1-p) u + (1-q) log(1-e^u)
	    */
	    w = (y == ML_NEGINF) // y = -Inf  well possible: we are on log scale!
		? 0. : (y - la) * exp(y - u + logbeta + r * u + t * R_Log1_Exp(u));
	    if(!R_FINITE(w))
		break;
	    if (i_pb >= n_N && w * wprev <= 0.)
		prev = fmax2(fabs(adj),fpu);
	    R_ifDEBUG_printf(
		"N(i=%2d): u=%#20.16g, pb(e^u)=%#15.9g, w=%#15.9g, %s prev=%g,",
		i_pb, u, y, w,
		(i_pb >= n_N && w * wprev <= 0.) ? "new" : "old", prev);
	    g = 1;
	    for (i_inn=0; i_inn < 1000; i_inn++) {
		adj = g * w;
		// safe guard (here, from the very beginning)
		if (fabs(adj) < prev) {
		    u_n = u - adj; // u_{n+1} = u_n - g*w
		    if (u_n <= 0.) { // <==> 0 <  xinbta := e^u  <= 1
			if (prev <= acu || fabs(w) <= acu) {
		 	    R_ifDEBUG_printf(
				" it{in}=%d, -adj=%g, %s <= acu  ==> convergence\n",
				i_inn, -adj, (prev <= acu) ? "prev" : "|w|");
			    goto L_converged;
			}
			// if (u_n != ML_NEGINF && u_n != 1)
			break;
		    }
		}
		g /= 3;
	    }
	    // (cancellation in (u_n -u) => may differ from adj:
	    double D = fmin2(fabs(adj), fabs(u_n - u));
	    /* R_ifDEBUG_printf(" delta(u)=%g\n", u_n - u); */
	    R_ifDEBUG_printf(" it{in}=%d, delta(u)=%9.3g, D/|.|=%.3g\n",
			     i_inn, u_n - u, D/fabs(u_n + u));
	    if (D <= 4e-16 * fabs(u_n + u))
		goto L_converged;
	    u = u_n;
	    xinbta = exp(u);
	    wprev = w;
	} // for(i )

    } else { // "normal scale" Newton

	for (i_pb=0; i_pb < 1000; i_pb++) {
	    y = pbeta_raw(xinbta, pp, qq, /*lower_tail = */ TRUE, log_p);
	    // delta{y} :   d_y = y - (log_p ? la : a);
#ifdef IEEE_754
	    if(!R_FINITE(y) && !(log_p && y == ML_NEGINF))// y = -Inf  is ok if(log_p)
#else
	    if (errno)
#endif
		{ // ML_ERR_return_NAN :
		    ML_ERROR(ME_DOMAIN, "");
		    qb[0] = qb[1] = ML_NAN; return;
		}


	    /* w := Newton step size  (F(.) - a) / F'(.)  or,
	     * --   log: (lF - la) / (F' / F) = exp(lF) * (lF - la) / F'
	     */
	    w = log_p
		? (y - la) * exp(y + logbeta + r * log(xinbta) + t * log1p(-xinbta))
		: (y - a)  * exp(    logbeta + r * log(xinbta) + t * log1p(-xinbta));
	    if (i_pb >= n_N && w * wprev <= 0.)
		prev = fmax2(fabs(adj),fpu);
	    R_ifDEBUG_printf(
		"N(i=%2d): x0=%#17.15g, pb(x0)=%#15.9g, w=%#15.9g, %s prev=%g,",
		i_pb, xinbta, y, w,
		(i_pb >= n_N && w * wprev <= 0.) ? "new" : "old", prev);
	    g = 1;
	    for (i_inn=0; i_inn < 1000;i_inn++) {
		adj = g * w;
		// take full Newton steps at the beginning; only then safe guard:
		if (i_pb < n_N || fabs(adj) < prev) {
		    tx = xinbta - adj; // x_{n+1} = x_n - g*w
		    if (0. <= tx && tx <= 1.) {
			if (prev <= acu || fabs(w) <= acu) {
			    R_ifDEBUG_printf(" it{in}=%d, delta(x)=%g, %s <= acu  ==> convergence\n",
					     i_inn, -adj, (prev <= acu) ? "prev" : "|w|");
			    goto L_converged;
			}
			if (tx != 0. && tx != 1)
			    break;
		    }
		}
		g /= 3;
	    }
	    R_ifDEBUG_printf(" it{in}=%d, delta(x)=%g\n", i_inn, tx - xinbta);
	    if (fabs(tx - xinbta) <= 4e-16 * (tx + xinbta)) // "<=" : (.) == 0
		goto L_converged;
	    xinbta = tx;
	    if(tx == 0) // "we have lost"
		break;
	    wprev = w;
	} // for( i_pb ..)

    } // end{else : normal scale Newton}

    /*-- NOT converged: Iteration count --*/
    warned = TRUE;
    ML_ERROR(ME_PRECISION, "qbeta");

L_converged:
    log_ = log_p || use_log_x; // only for printing
    R_ifDEBUG_printf(" %s: Final delta(y) = %g%s\n",
		     warned ? "_NO_ convergence" : "converged",
		     y - (log_ ? la : a), (log_ ? " (log_)" : ""));
    if((log_ && y == ML_NEGINF) || (!log_ && y == 0)) {
	// stuck at left, try if smallest positive number is "better"
	w = pbeta_raw(DBL_very_MIN, pp, qq, TRUE, log_);
	if(log_ || fabs(w - a) <= fabs(y - a)) {
	    tx  = DBL_very_MIN;
	    u_n = DBL_log_v_MIN;// = log(DBL_very_MIN)
	}
	add_N_step = FALSE; // not trying to do better anymore
    }
    else if(!warned && (log_ ? fabs(y - la) > 3 : fabs(y - a) > 1e-4)) {
	if(!(log_ && y == ML_NEGINF &&
	     // e.g. qbeta(-1e-10, .2, .03, log=TRUE) cannot get accurate ==> do NOT warn
	     pbeta_raw(DBL_1__eps, // = 1 - eps
		       pp, qq, TRUE, TRUE) > la + 2))
	    MATHLIB_WARNING2( // low accuracy for more platform independent output:
		"qbeta(a, *) =: x0 with |pbeta(x0,*%s) - alpha| = %.5g is not accurate",
		(log_ ? ", log_" : ""), fabs(y - (log_ ? la : a)));
    }
L_return:
    if(give_log_q) { // ==> use_log_x , too
	if(!use_log_x) // (see if claim above is true)
	    MATHLIB_WARNING(
		"qbeta() L_return, u_n=%g;  give_log_q=TRUE but use_log_x=FALSE -- please report!",
		u_n);
	double r = R_Log1_Exp(u_n);
	if(swap_tail) {
	    qb[0] = r;	 qb[1] = u_n;
	} else {
	    qb[0] = u_n; qb[1] = r;
	}
    } else {
	if(use_log_x) {
	    if(add_N_step) {
		/* add one last Newton step on original x scale, e.g., for
		   qbeta(2^-98, 0.125, 2^-96) */
		xinbta = exp(u_n);
		y = pbeta_raw(xinbta, pp, qq, /*lower_tail = */ TRUE, log_p);
		w = log_p
		    ? (y - la) * exp(y + logbeta + r * log(xinbta) + t * log1p(-xinbta))
		    : (y - a)  * exp(    logbeta + r * log(xinbta) + t * log1p(-xinbta));
		tx = xinbta - w;
		R_ifDEBUG_printf(" Final Newton correction(non-log scale):\n"
								   //   \n  xinbta=%.16g
				 "  xinbta=%.16g, y=%g, w=-Delta(x)=%g. \n=> new x=%.16g\n",
		    xinbta, y, w, tx);
	    } else {
		if(swap_tail) {
		    qb[0] = -expm1(u_n); qb[1] =  exp  (u_n);
		} else {
		    qb[0] =  exp  (u_n); qb[1] = -expm1(u_n);
		}
		return;
	    }
	}
	if(swap_tail) {
	    qb[0] = 1 - tx; qb[1] = tx;
	} else {
	    qb[0] = tx;	qb[1] = 1 - tx;
	}
    }
    return;
}
/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 2006 The R Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 *
 *  SYNOPSIS
 *
 * #include <Rmath.h>
 *
 * double pbeta_raw(double x, double a, double b, int lower_tail, int log_p)
 * double pbeta	   (double x, double a, double b, int lower_tail, int log_p)
 *
 *  DESCRIPTION
 *
 *	Returns distribution function of the beta distribution.
 *	( = The incomplete beta ratio I_x(p,q) ).
 *
 *  NOTES
 *
 *      As from R 2.3.0, a wrapper for TOMS708
 *      as from R 2.6.0, 'log_p' partially improved over log(p..)
 */

#include "from-R.h"


attribute_hidden
double pbeta_raw(double x, double a, double b, int lower_tail, int log_p)
{
    // treat limit cases correctly here:
    if(a == 0 || b == 0 || !R_FINITE(a) || !R_FINITE(b)) {
	// NB:  0 < x < 1 :
	if(a == 0 && b == 0) // point mass 1/2 at each of {0,1} :
	    return (log_p ? -M_LN2 : 0.5);
	if (a == 0 || a/b == 0) // point mass 1 at 0 ==> P(X <= x) = 1, all x > 0
	    return R_DT_1;
	if (b == 0 || b/a == 0) // point mass 1 at 1 ==> P(X <= x) = 0, all x < 1
	    return R_DT_0;
	// else, remaining case:  a = b = Inf : point mass 1 at 1/2
	if (x < 0.5) return R_DT_0; else return R_DT_1;
    }
    // Now:  0 < a < Inf;  0 < b < Inf

    double x1 = 0.5 - x + 0.5, w, wc;
    int ierr;
    //====
    bratio(a, b, x, x1, &w, &wc, &ierr, log_p); /* -> ./toms708.c */
    //====
    // ierr in {10,14} <==> bgrat() error code ierr-10 in 1:4; for 1 and 4, warned *there*
    if(ierr && ierr != 11 && ierr != 14)
	MATHLIB_WARNING4(_("pbeta_raw(%g, a=%g, b=%g, ..) -> bratio() gave error code %d"),
			x, a,b, ierr);
    return lower_tail ? w : wc;
} /* pbeta_raw() */

double pbeta(double x, double a, double b, int lower_tail, int log_p)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(a) || ISNAN(b)) return x + a + b;
#endif

    if (a < 0 || b < 0) ML_ERR_return_NAN;
    // allowing a==0 and b==0  <==> treat as one- or two-point mass

    if (x <= 0)
	return R_DT_0;
    if (x >= 1)
	return R_DT_1;

    return pbeta_raw(x, a, b, lower_tail, log_p);
}
// -*- mode: C ; delete-old-versions: never -*-

/* Based on C translation of ACM TOMS 708
   Please do not change this, e.g. to use R's versions of the
   ancillary routines, without investigating the error analysis as we
   do need very high relative accuracy.  This version has about
   14 digits accuracy.
*/

#undef min
#define min(a,b) ((a < b)?a:b)
#undef max
#define max(a,b) ((a > b)?a:b)

#include "from-R.h" /* includes config.h, math.h */

/* after config.h to avoid warning on Solaris */
#include <limits.h>

/**----------- DEBUGGING -------------
 *
 *	make CFLAGS='-DDEBUG_bratio  ...'
 *MM (w/ Debug, w/o Optimization):
 (cd `R-devel-pbeta-dbg RHOME`/src/nmath ; gcc -I. -I../../src/include -I../../../R/src/include  -DHAVE_CONFIG_H -fopenmp -g -pedantic -Wall --std=gnu99 -DDEBUG_q -DDEBUG_bratio -Wcast-align -Wclobbered  -c ../../../R/src/nmath/toms708.c -o toms708.o; cd ../..; make R)
*/
#ifdef DEBUG_bratio
# define R_ifDEBUG_printf(...) REprintf(__VA_ARGS__)
#else
# define R_ifDEBUG_printf(...)
#endif

/* MM added R_D_LExp, so redefine here in terms of rexpm1 */
#undef R_Log1_Exp
#define R_Log1_Exp(x)   ((x) > -M_LN2 ? log(-rexpm1(x)) : log1p(-exp(x)))


static double bfrac(double, double, double, double, double, double, int log_p);
static void bgrat(double, double, double, double, double *, double, int *, Rboolean log_w);
static double grat_r(double a, double x, double r, double eps);
static double apser(double, double, double, double);
static double bpser(double, double, double, double, int log_p);
static double basym(double, double, double, double, int log_p);
static double fpser(double, double, double, double, int log_p);
static double bup(double, double, double, double, int, double, int give_log);
static double exparg(int);
static double psi(double);
static double gam1(double);
static double gamln1(double);
static double betaln(double, double);
static double algdiv(double, double);
static double brcmp1(int, double, double, double, double, int give_log);
static double brcomp(double, double, double, double, int log_p);
static double rlog1(double);
static double bcorr(double, double);
static double gamln(double);
static double alnrel(double);
static double esum(int, double, int give_log);
static double erf__(double);
static double rexpm1(double);
static double erfc1(int, double);
static double gsumln(double, double);

/*      ALGORITHM 708, COLLECTED ALGORITHMS FROM ACM.
 *      This work published in  Transactions On Mathematical Software,
 *      vol. 18, no. 3, September 1992, pp. 360-373.
 */

/* Changes by R Core Team :
 * add log_p  and work towards gaining precision in that case
 */

void attribute_hidden
bratio(double a, double b, double x, double y, double *w, double *w1,
       int *ierr, int log_p)
{
/* -----------------------------------------------------------------------

 *	      Evaluation of the Incomplete Beta function I_x(a,b)

 *		       --------------------

 *     It is assumed that a and b are nonnegative, and that x <= 1
 *     and y = 1 - x.  Bratio assigns w and w1 the values

 *			w  = I_x(a,b)
 *			w1 = 1 - I_x(a,b)

 *     ierr is a variable that reports the status of the results.
 *     If no input errors are detected then ierr is set to 0 and
 *     w and w1 are computed. otherwise, if an error is detected,
 *     then w and w1 are assigned the value 0 and ierr is set to
 *     one of the following values ...

 *	  ierr = 1  if a or b is negative
 *	  ierr = 2  if a = b = 0
 *	  ierr = 3  if x < 0 or x > 1
 *	  ierr = 4  if y < 0 or y > 1
 *	  ierr = 5  if x + y != 1
 *	  ierr = 6  if x = a = 0
 *	  ierr = 7  if y = b = 0
 *	  ierr = 8	(not used currently)
 *	  ierr = 9  NaN in a, b, x, or y
 *	  ierr = 10     (not used currently)
 *	  ierr = 11  bgrat() error code 1 [+ warning in bgrat()]
 *	  ierr = 12  bgrat() error code 2   (no warning here)
 *	  ierr = 13  bgrat() error code 3   (no warning here)
 *	  ierr = 14  bgrat() error code 4 [+ WARNING in bgrat()]


 * --------------------
 *     Written by Alfred H. Morris, Jr.
 *	  Naval Surface Warfare Center
 *	  Dahlgren, Virginia
 *     Revised ... Nov 1991
* ----------------------------------------------------------------------- */

    Rboolean do_swap;
    int n, ierr1 = 0;
    double z, a0, b0, x0, y0, lambda;

/*  eps is a machine dependent constant: the smallest
 *      floating point number for which   1. + eps > 1.
 * NOTE: for almost all purposes it is replaced by 1e-15 (~= 4.5 times larger) below */
    double eps = 2. * Rf_d1mach(3); /* == DBL_EPSILON (in R, Rmath) */

/* ----------------------------------------------------------------------- */
    *w  = R_D__0;
    *w1 = R_D__0;

#ifdef IEEE_754
    // safeguard, preventing infinite loops further down
    if (ISNAN(x) || ISNAN(y) ||
	ISNAN(a) || ISNAN(b)) { *ierr = 9; return; }
#endif
    if (a < 0. || b < 0.)   { *ierr = 1; return; }
    if (a == 0. && b == 0.) { *ierr = 2; return; }
    if (x < 0. || x > 1.)   { *ierr = 3; return; }
    if (y < 0. || y > 1.)   { *ierr = 4; return; }

    /* check that  'y == 1 - x' : */
    z = x + y - 0.5 - 0.5;

    if (fabs(z) > eps * 3.) { *ierr = 5; return; }

    R_ifDEBUG_printf("bratio(a=%g, b=%g, x=%9g, y=%9g, .., log_p=%d): ",
		     a,b,x,y, log_p);
    *ierr = 0;
    if (x == 0.) goto L200;
    if (y == 0.) goto L210;

    if (a == 0.) goto L211;
    if (b == 0.) goto L201;

    eps = max(eps, 1e-15);
    Rboolean a_lt_b = (a < b);
    if (/* max(a,b) */ (a_lt_b ? b : a) < eps * .001) { /* procedure for a and b < 0.001 * eps */
	// L230:  -- result *independent* of x (!)
	// *w  = a/(a+b)  and  w1 = b/(a+b) :
	if(log_p) {
	    if(a_lt_b) {
		*w  = log1p(-a/(a+b)); // notably if a << b
		*w1 = log  ( a/(a+b));
	    } else { // b <= a
		*w  = log  ( b/(a+b));
		*w1 = log1p(-b/(a+b));
	    }
	} else {
	    *w	= b / (a + b);
	    *w1 = a / (a + b);
	}

	R_ifDEBUG_printf("a & b very small -> simple ratios (%g,%g)\n", *w,*w1);
	return;
    }

#define SET_0_noswap \
    a0 = a;  x0 = x; \
    b0 = b;  y0 = y;

#define SET_0_swap   \
    a0 = b;  x0 = y; \
    b0 = a;  y0 = x;

    if (min(a,b) <= 1.) { /*------------------------ a <= 1  or  b <= 1 ---- */

	do_swap = (x > 0.5);
	if (do_swap) {
	    SET_0_swap;
	} else {
	    SET_0_noswap;
	}
	/* now have  x0 <= 1/2 <= y0  (still  x0+y0 == 1) */

	R_ifDEBUG_printf(" min(a,b) <= 1, do_swap=%d;", do_swap);

	if (b0 < min(eps, eps * a0)) { /* L80: */
	    *w = fpser(a0, b0, x0, eps, log_p);
	    *w1 = log_p ? R_Log1_Exp(*w) : 0.5 - *w + 0.5;
	    R_ifDEBUG_printf("  b0 small -> w := fpser(*) = %.15g\n", *w);
	    goto L_end;
	}

	if (a0 < min(eps, eps * b0) && b0 * x0 <= 1.) { /* L90: */
	    *w1 = apser(a0, b0, x0, eps);
	    R_ifDEBUG_printf("  a0 small -> w1 := apser(*) = %.15g\n", *w1);
	    goto L_end_from_w1;
	}

	Rboolean did_bup = FALSE;
	if (max(a0,b0) > 1.) { /* L20:  min(a,b) <= 1 < max(a,b)  */
	    R_ifDEBUG_printf("\n L20:  min(a,b) <= 1 < max(a,b); ");
	    if (b0 <= 1.) goto L_w_bpser;

	    if (x0 >= 0.29) /* was 0.3, PR#13786 */	goto L_w1_bpser;

	    if (x0 < 0.1 && pow(x0*b0, a0) <= 0.7)	goto L_w_bpser;

	    if (b0 > 15.) {
		*w1 = 0.;
		goto L131;
	    }
	} else { /*  a, b <= 1 */
	    R_ifDEBUG_printf("\n      both a,b <= 1; ");
	    if (a0 >= min(0.2, b0))	goto L_w_bpser;

	    if (pow(x0, a0) <= 0.9) 	goto L_w_bpser;

	    if (x0 >= 0.3)		goto L_w1_bpser;
	}
	n = 20; /* goto L130; */
	*w1 = bup(b0, a0, y0, x0, n, eps, FALSE); did_bup = TRUE;
	R_ifDEBUG_printf("  ... n=20 and *w1 := bup(*) = %.15g; ", *w1);
	b0 += n;
    L131:
	R_ifDEBUG_printf(" L131: bgrat(*, w1=%.15g) ", *w1);
	bgrat(b0, a0, y0, x0, w1, 15*eps, &ierr1, FALSE);
#ifdef DEBUG_bratio
	REprintf(" ==> new w1=%.15g", *w1);
	if(ierr1) REprintf(" ERROR(code=%d)\n", ierr1) ; else REprintf("\n");
#endif
	if(*w1 == 0 || (0 < *w1 && *w1 < DBL_MIN)) { // w1=0 or very close:
	    // "almost surely" from underflow, try more: [2013-03-04]
// FIXME: it is even better to do this in bgrat *directly* at least for the case
//  !did_bup, i.e., where *w1 = (0 or -Inf) on entry
	    R_ifDEBUG_printf(" denormalized or underflow (?) -> retrying: ");
	    if(did_bup) { // re-do that part on log scale:
		*w1 = bup(b0-n, a0, y0, x0, n, eps, TRUE);
	    }
	    else *w1 = ML_NEGINF; // = 0 on log-scale
	    bgrat(b0, a0, y0, x0, w1, 15*eps, &ierr1, TRUE);
	    if(ierr1) *ierr = 10 + ierr1;
#ifdef DEBUG_bratio
	    REprintf(" ==> new log(w1)=%.15g", *w1);
	    if(ierr1) REprintf(" Error(code=%d)\n", ierr1) ; else REprintf("\n");
#endif
	    goto L_end_from_w1_log;
	}
	// else
	if(ierr1) *ierr = 10 + ierr1;
	if(*w1 < 0)
	    MATHLIB_WARNING4("bratio(a=%g, b=%g, x=%g): bgrat() -> w1 = %g",
			     a,b,x, *w1);
	goto L_end_from_w1;
    }
    else { /* L30: -------------------- both  a, b > 1  {a0 > 1  &  b0 > 1} ---*/

	/* lambda := a y - b x  =  (a + b)y  =  a - (a+b)x    {using x + y == 1},
	 * ------ using the numerically best version : */
	lambda = R_FINITE(a+b)
	    ? ((a > b) ? (a + b) * y - b : a - (a + b) * x)
	    : a*y - b*x;
	do_swap = (lambda < 0.);
	if (do_swap) {
	    lambda = -lambda;
	    SET_0_swap;
	} else {
	    SET_0_noswap;
	}

	R_ifDEBUG_printf("  L30:  both  a, b > 1; |lambda| = %#g, do_swap = %d\n",
			 lambda, do_swap);

	if (b0 < 40.) {
	    R_ifDEBUG_printf("  b0 < 40;");
	    if (b0 * x0 <= 0.7
		|| (log_p && lambda > 650.)) // << added 2010-03; svn r51327
		goto L_w_bpser;
	    else
		goto L140;
	}
	else if (a0 > b0) { /* ----  a0 > b0 >= 40  ---- */
	    R_ifDEBUG_printf("  a0 > b0 >= 40;");
	    if (b0 <= 100. || lambda > b0 * 0.03)
		goto L_bfrac;

	} else if (a0 <= 100.) {
	    R_ifDEBUG_printf("  a0 <= 100; a0 <= b0 >= 40;");
	    goto L_bfrac;
	}
	else if (lambda > a0 * 0.03) {
	    R_ifDEBUG_printf("  b0 >= a0 > 100; lambda > a0 * 0.03 ");
	    goto L_bfrac;
	}

	/* else if none of the above    L180: */
	*w = basym(a0, b0, lambda, eps * 100., log_p);
	*w1 = log_p ? R_Log1_Exp(*w) : 0.5 - *w + 0.5;
	R_ifDEBUG_printf("  b0 >= a0 > 100; lambda <= a0 * 0.03: *w:= basym(*) =%.15g\n",
			 *w);
	goto L_end;

    } /* else: a, b > 1 */

/*            EVALUATION OF THE APPROPRIATE ALGORITHM */

L_w_bpser: // was L100
    *w = bpser(a0, b0, x0, eps, log_p);
    *w1 = log_p ? R_Log1_Exp(*w) : 0.5 - *w + 0.5;
    R_ifDEBUG_printf(" L_w_bpser: *w := bpser(*) = %.15g\n", *w);
    goto L_end;

L_w1_bpser:  // was L110
    *w1 = bpser(b0, a0, y0, eps, log_p);
    *w  = log_p ? R_Log1_Exp(*w1) : 0.5 - *w1 + 0.5;
    R_ifDEBUG_printf(" L_w1_bpser: *w1 := bpser(*) = %.15g\n", *w1);
    goto L_end;

L_bfrac:
    *w = bfrac(a0, b0, x0, y0, lambda, eps * 15., log_p);
    *w1 = log_p ? R_Log1_Exp(*w) : 0.5 - *w + 0.5;
    R_ifDEBUG_printf(" L_bfrac: *w := bfrac(*) = %g\n", *w);
    goto L_end;

L140:
    /* b0 := fractional_part( b0 )  in (0, 1]  */
    n = (int) b0;
    b0 -= n;
    if (b0 == 0.) {
	--n; b0 = 1.;
    }

    *w = bup(b0, a0, y0, x0, n, eps, FALSE);

    if(*w < DBL_MIN && log_p) { /* do not believe it; try bpser() : */
	R_ifDEBUG_printf(" L140: bup(b0=%g,..)=%.15g < DBL_MIN - not used; ", b0, *w);
	/*revert: */ b0 += n;
	/* which is only valid if b0 <= 1 || b0*x0 <= 0.7 */
	goto L_w_bpser;
    }
    R_ifDEBUG_printf(" L140: *w := bup(b0=%g,..) = %.15g; ", b0, *w);
    if (x0 <= 0.7) {
	/* log_p :  TODO:  w = bup(.) + bpser(.)  -- not so easy to use log-scale */
	*w += bpser(a0, b0, x0, eps, /* log_p = */ FALSE);
	R_ifDEBUG_printf(" x0 <= 0.7: *w := *w + bpser(*) = %.15g\n", *w);
	goto L_end_from_w;
    }
    /* L150: */
    if (a0 <= 15.) {
	n = 20;
	*w += bup(a0, b0, x0, y0, n, eps, FALSE);
	R_ifDEBUG_printf("\n a0 <= 15: *w := *w + bup(*) = %.15g;", *w);
	a0 += n;
    }
    R_ifDEBUG_printf(" bgrat(*, w=%.15g) ", *w);
    bgrat(a0, b0, x0, y0, w, 15*eps, &ierr1, FALSE);
    if(ierr1) *ierr = 10 + ierr1;
#ifdef DEBUG_bratio
    REprintf("==> new w=%.15g", *w);
    if(ierr1) REprintf(" Error(code=%d)\n", ierr1) ; else REprintf("\n");
#endif
    goto L_end_from_w;


/* TERMINATION OF THE PROCEDURE */

L200:
    if (a == 0.) { *ierr = 6;    return; }
    // else:
L201: *w  = R_D__0; *w1 = R_D__1; return;

L210:
    if (b == 0.) { *ierr = 7;    return; }
    // else:
L211: *w  = R_D__1; *w1 = R_D__0; return;

L_end_from_w:
    if(log_p) {
	*w1 = log1p(-*w);
	*w  = log(*w);
    } else {
	*w1 = 0.5 - *w + 0.5;
    }
    goto L_end;

L_end_from_w1:
    if(log_p) {
	*w  = log1p(-*w1);
	*w1 = log(*w1);
    } else {
	*w = 0.5 - *w1 + 0.5;
    }
    goto L_end;

L_end_from_w1_log:
    // *w1 = log(w1) already; w = 1 - w1  ==> log(w) = log(1 - w1) = log(1 - exp(*w1))
    if(log_p) {
	*w = R_Log1_Exp(*w1);
    } else {
	*w  = /* 1 - exp(*w1) */ -expm1(*w1);
	*w1 = exp(*w1);
    }
    goto L_end;


L_end:
    if (do_swap) { /* swap */
	double t = *w; *w = *w1; *w1 = t;
    }
    return;

} /* bratio */

#undef SET_0_noswap
#undef SET_0_swap

double fpser(double a, double b, double x, double eps, int log_p)
{
/* ----------------------------------------------------------------------- *

 *                 EVALUATION OF I (A,B)
 *                                X

 *          FOR B < MIN(EPS, EPS*A) AND X <= 0.5

 * ----------------------------------------------------------------------- */

    double ans, c, s, t, an, tol;

    /* SET  ans := x^a : */
    if (log_p) {
	ans = a * log(x);
    } else if (a > eps * 0.001) {
	t = a * log(x);
	if (t < exparg(1)) { /* exp(t) would underflow */
	    return 0.;
	}
	ans = exp(t);
    } else
	ans = 1.;

/*                NOTE THAT 1/B(A,B) = B */

    if (log_p)
	ans += log(b) - log(a);
    else
	ans *= b / a;

    tol = eps / a;
    an = a + 1.;
    t = x;
    s = t / an;
    do {
	an += 1.;
	t = x * t;
	c = t / an;
	s += c;
    } while (fabs(c) > tol);

    if (log_p)
	ans += log1p(a * s);
    else
	ans *= a * s + 1.;
    return ans;
} /* fpser */

static double apser(double a, double b, double x, double eps)
{
/* -----------------------------------------------------------------------
 *     apser() yields the incomplete beta ratio  I_{1-x}(b,a)  for
 *     a <= min(eps,eps*b), b*x <= 1, and x <= 0.5,  i.e., a is very small.
 *     Use only if above inequalities are satisfied.
 * ----------------------------------------------------------------------- */

    static double const g = .577215664901533;

    double tol, c, j, s, t, aj;
    double bx = b * x;

    t = x - bx;
    if (b * eps <= 0.02)
	c = log(x) + psi(b) + g + t;
    else // b > 2e13 : psi(b) ~= log(b)
	c = log(bx) + g + t;

    tol = eps * 5. * fabs(c);
    j = 1.;
    s = 0.;
    do {
	j += 1.;
	t *= x - bx / j;
	aj = t / j;
	s += aj;
    } while (fabs(aj) > tol);

    return -a * (c + s);
} /* apser */

static double bpser(double a, double b, double x, double eps, int log_p)
{
/* -----------------------------------------------------------------------
 * Power SERies expansion for evaluating I_x(a,b) when
 *	       b <= 1 or b*x <= 0.7.   eps is the tolerance used.
 * NB: if log_p is TRUE, also use it if   (b < 40  & lambda > 650)
 * ----------------------------------------------------------------------- */

    int i, m;
    double ans, c, t, u, z, a0, b0, apb;

    if (x == 0.) {
	return R_D__0;
    }
/* ----------------------------------------------------------------------- */
/*	      compute the factor  x^a/(a*Beta(a,b)) */
/* ----------------------------------------------------------------------- */
    a0 = min(a,b);
    if (a0 >= 1.) { /*		 ------	 1 <= a0 <= b0  ------ */
	z = a * log(x) - betaln(a, b);
	ans = log_p ? z - log(a) : exp(z) / a;
    }
    else {
	b0 = max(a,b);

	if (b0 < 8.) {

	    if (b0 <= 1.) { /*	 ------	 a0 < 1	 and  b0 <= 1  ------ */

		if(log_p) {
		    ans = a * log(x);
		} else {
		    ans = pow(x, a);
		    if (ans == 0.) /* once underflow, always underflow .. */
			return ans;
		}
		apb = a + b;
		if (apb > 1.) {
		    u = a + b - 1.;
		    z = (gam1(u) + 1.) / apb;
		} else {
		    z = gam1(apb) + 1.;
		}
		c = (gam1(a) + 1.) * (gam1(b) + 1.) / z;

		if(log_p) /* FIXME ? -- improve quite a bit for c ~= 1 */
		    ans += log(c * (b / apb));
		else
		    ans *=  c * (b / apb);

	    } else { /* 	------	a0 < 1 < b0 < 8	 ------ */

		u = gamln1(a0);
		m = (int)(b0 - 1.);
		if (m >= 1) {
		    c = 1.;
		    for (i = 1; i <= m; ++i) {
			b0 += -1.;
			c *= b0 / (a0 + b0);
		    }
		    u += log(c);
		}

		z = a * log(x) - u;
		b0 += -1.; // => b0 in (0, 7)
		apb = a0 + b0;
		if (apb > 1.) {
		    u = a0 + b0 - 1.;
		    t = (gam1(u) + 1.) / apb;
		} else {
		    t = gam1(apb) + 1.;
		}

		if(log_p) /* FIXME? potential for improving log(t) */
		    ans = z + log(a0 / a) + log1p(gam1(b0)) - log(t);
		else
		    ans = exp(z) * (a0 / a) * (gam1(b0) + 1.) / t;
	    }

	} else { /* 		------  a0 < 1 < 8 <= b0  ------ */

	    u = gamln1(a0) + algdiv(a0, b0);
	    z = a * log(x) - u;

	    if(log_p)
		ans = z + log(a0 / a);
	    else
		ans = a0 / a * exp(z);
	}
    }
    R_ifDEBUG_printf(" bpser(a=%g, b=%g, x=%g, log=%d): prelim.ans = %.14g;\n",
		     a,b,x, log_p, ans);
    if (ans == R_D__0 || (!log_p && a <= eps * 0.1)) {
	return ans;
    }

/* ----------------------------------------------------------------------- */
/*		       COMPUTE THE SERIES */
/* ----------------------------------------------------------------------- */
    double tol = eps / a,
	n = 0.,
	sum = 0., w;
    c = 1.;
    do { // sum is alternating as long as n < b (<==> 1 - b/n < 0)
	n += 1.;
	c *= (0.5 - b / n + 0.5) * x;
	w = c / (a + n);
	sum += w;
    } while (n < 1e7 && fabs(w) > tol);
    if(fabs(w) > tol) { // the series did not converge (in time)
	// warn only when the result seems to matter:
	if(( log_p && !(a*sum > -1. && fabs(log1p(a * sum)) < eps*fabs(ans))) ||
	   (!log_p && fabs(a*sum + 1.) != 1.))
	    MATHLIB_WARNING5(
		" bpser(a=%g, b=%g, x=%g,...) did not converge (n=1e7, |w|/tol=%g > 1; A=%g)",
		a,b,x, fabs(w)/tol, ans);
    }
    R_ifDEBUG_printf("  -> n=%.0f iterations, |w|=%g %s %g=tol:=eps/a ==> a*sum=%g\n",
		     n, fabs(w), (fabs(w) > tol) ? ">!!>" : "<=",
		     tol, a*sum);
    if(log_p) {
	if (a*sum > -1.) ans += log1p(a * sum);
	else {
	    if(ans > ML_NEGINF)
		MATHLIB_WARNING3(
		    "pbeta(*, log.p=TRUE) -> bpser(a=%g, b=%g, x=%g,...) underflow to -Inf",
		    a,b,x);
	    ans = ML_NEGINF;
	}
    } else if (a*sum > -1.)
	ans *= (a * sum + 1.);
    else // underflow to
	ans = 0.;
    return ans;
} /* bpser */

static double bup(double a, double b, double x, double y, int n, double eps,
		  int give_log)
{
/* ----------------------------------------------------------------------- */
/*     EVALUATION OF I_x(A,B) - I_x(A+N,B) WHERE N IS A POSITIVE INT. */
/*     EPS IS THE TOLERANCE USED. */
/* ----------------------------------------------------------------------- */

    double ret_val;
    int i, k, mu;
    double d, l;

// Obtain the scaling factor exp(-mu) and exp(mu)*(x^a * y^b / beta(a,b))/a

    double apb = a + b,
	ap1 = a + 1.;
    if (n > 1 && a >= 1. && apb >= ap1 * 1.1) {
	mu = (int)fabs(exparg(1));
	k = (int) exparg(0);
	if (mu > k)
	    mu = k;
 	d = exp(-(double) mu);
    }
    else {
	mu = 0;
	d = 1.;
    }

    /* L10: */
    ret_val = give_log
	? brcmp1(mu, a, b, x, y, TRUE) - log(a)
	: brcmp1(mu, a, b, x, y, FALSE)  / a;
    if (n == 1 ||
	(give_log && ret_val == ML_NEGINF) || (!give_log && ret_val == 0.))
	return ret_val;

    int nm1 = n - 1;
    double w = d;

/*          LET K BE THE INDEX OF THE MAXIMUM TERM */

    k = 0;
    if (b > 1.) {
	if (y > 1e-4) {
	    double r = (b - 1.) * x / y - a;
	    if (r >= 1.)
		k = (r < nm1) ? (int) r : nm1;
	} else
	    k = nm1;

//          ADD THE INCREASING TERMS OF THE SERIES - if k > 0
/* L30: */
	for (i = 0; i < k; ++i) {
	    l = (double) i;
	    d *= (apb + l) / (ap1 + l) * x;
	    w += d;
	}
    }

// L40:     ADD THE REMAINING TERMS OF THE SERIES

    for (i = k; i < nm1; ++i) {
	l = (double) i;
	d *= (apb + l) / (ap1 + l) * x;
	w += d;
	if (d <= eps * w) /* relativ convergence (eps) */
	    break;
    }

    // L50: TERMINATE THE PROCEDURE
    if(give_log) {
	ret_val += log(w);
    } else
	ret_val *= w;

    return ret_val;
} /* bup */

static double bfrac(double a, double b, double x, double y, double lambda,
		    double eps, int log_p)
{
/* -----------------------------------------------------------------------
       Continued fraction expansion for I_x(a,b) when a, b > 1.
       It is assumed that  lambda = (a + b)*y - b.
   -----------------------------------------------------------------------*/

    double c, e, n, p, r, s, t, w, c0, c1, r0, an, bn, yp1, anp1, bnp1,
	beta, alpha, brc;

    if(!R_FINITE(lambda)) return ML_NAN;// TODO: can return 0 or 1 (?)
    R_ifDEBUG_printf(" bfrac(a=%g, b=%g, x=%g, y=%g, lambda=%g, eps=%g, log_p=%d):",
		     a,b,x,y, lambda, eps, log_p);
    brc = brcomp(a, b, x, y, log_p);
    if(ISNAN(brc)) { // e.g. from   L <- 1e308; pnbinom(L, L, mu = 5)
	R_ifDEBUG_printf(" --> brcomp(a,b,x,y) = NaN\n");
	ML_ERR_return_NAN; // TODO: could we know better?
    }
    if (!log_p && brc == 0.) {
	R_ifDEBUG_printf(" --> brcomp(a,b,x,y) underflowed to 0.\n");
	return 0.;
    }
#ifdef DEBUG_bratio
    else
	REprintf("\n");
#endif

    c = lambda + 1.;
    c0 = b / a;
    c1 = 1. / a + 1.;
    yp1 = y + 1.;

    n = 0.;
    p = 1.;
    s = a + 1.;
    an = 0.;
    bn = 1.;
    anp1 = 1.;
    bnp1 = c / c1;
    r = c1 / c;

/*        CONTINUED FRACTION CALCULATION */

    do {
	n += 1.;
	t = n / a;
	w = n * (b - n) * x;
	e = a / s;
	alpha = p * (p + c0) * e * e * (w * x);
	e = (t + 1.) / (c1 + t + t);
	beta = n + w / s + e * (c + n * yp1);
	p = t + 1.;
	s += 2.;

	/* update an, bn, anp1, and bnp1 */

	t = alpha * an + beta * anp1;	an = anp1;	anp1 = t;
	t = alpha * bn + beta * bnp1;	bn = bnp1;	bnp1 = t;

	r0 = r;
	r = anp1 / bnp1;
#ifdef _not_normally_DEBUG_bfrac
	R_ifDEBUG_printf(" n=%5.0f, a_{n,n+1}= (%12g,%12g),  b_{n,n+1} = (%12g,%12g) => r0,r = (%14g,%14g)\n",
			 n, an,anp1, bn,bnp1, r0, r);
#endif
	if (fabs(r - r0) <= eps * r)
	    break;

	/* rescale an, bn, anp1, and bnp1 */

	an /= bnp1;
	bn /= bnp1;
	anp1 = r;
	bnp1 = 1.;
    } while (n < 10000);// arbitrary; had '1' --> infinite loop for  lambda = Inf
    R_ifDEBUG_printf("  in bfrac(): n=%.0f terms cont.frac.; brc=%g, r=%g\n",
		     n, brc, r);
    if(n >= 10000 && fabs(r - r0) > eps * r)
	MATHLIB_WARNING5(
	    " bfrac(a=%g, b=%g, x=%g, y=%g, lambda=%g) did *not* converge (in 10000 steps)\n",
	    a,b,x,y, lambda);
    return (log_p ? brc + log(r) : brc * r);
} /* bfrac */

static double brcomp(double a, double b, double x, double y, int log_p)
{
/* -----------------------------------------------------------------------
 *		 Evaluation of x^a * y^b / Beta(a,b)
 * ----------------------------------------------------------------------- */

    static double const__ = .398942280401433; /* == 1/sqrt(2*pi); */
    /* R has  M_1_SQRT_2PI , and M_LN_SQRT_2PI = ln(sqrt(2*pi)) = 0.918938.. */
    int i, n;
    double c, e, u, v, z, a0, b0, apb;

    if (x == 0. || y == 0.) {
	return R_D__0;
    }
    a0 = min(a, b);
    if (a0 < 8.) {
	double lnx, lny;
	if (x <= .375) {
	    lnx = log(x);
	    lny = alnrel(-x);
	}
	else {
	    if (y > .375) {
		lnx = log(x);
		lny = log(y);
	    } else {
		lnx = alnrel(-y);
		lny = log(y);
	    }
	}

	z = a * lnx + b * lny;
	if (a0 >= 1.) {
	    z -= betaln(a, b);
	    return R_D_exp(z);
	}

/* ----------------------------------------------------------------------- */
/*		PROCEDURE FOR a < 1 OR b < 1 */
/* ----------------------------------------------------------------------- */

	b0 = max(a, b);
	if (b0 >= 8.) { /* L80: */
	    u = gamln1(a0) + algdiv(a0, b0);

	    return (log_p ? log(a0) + (z - u)  : a0 * exp(z - u));
	}
	/* else : */

	if (b0 <= 1.) { /*		algorithm for max(a,b) = b0 <= 1 */

	    double e_z = R_D_exp(z);

	    if (!log_p && e_z == 0.) /* exp() underflow */
		return 0.;

	    apb = a + b;
	    if (apb > 1.) {
		u = a + b - 1.;
		z = (gam1(u) + 1.) / apb;
	    } else {
		z = gam1(apb) + 1.;
	    }

	    c = (gam1(a) + 1.) * (gam1(b) + 1.) / z;
	    /* FIXME? log(a0*c)= log(a0)+ log(c) and that is improvable */
	    return (log_p
		    ? e_z + log(a0 * c) - log1p(a0/b0)
		    : e_z * (a0 * c) / (a0 / b0 + 1.));
	}

	/* else : 		  ALGORITHM FOR 1 < b0 < 8 */

	u = gamln1(a0);
	n = (int)(b0 - 1.);
	if (n >= 1) {
	    c = 1.;
	    for (i = 1; i <= n; ++i) {
		b0 += -1.;
		c *= b0 / (a0 + b0);
	    }
	    u = log(c) + u;
	}
	z -= u;
	b0 += -1.;
	apb = a0 + b0;
	double t;
	if (apb > 1.) {
	    u = a0 + b0 - 1.;
	    t = (gam1(u) + 1.) / apb;
	} else {
	    t = gam1(apb) + 1.;
	}

	return (log_p
		? log(a0) + z + log1p(gam1(b0))  - log(t)
		: a0 * exp(z) * (gam1(b0) + 1.) / t);

    } else {
/* ----------------------------------------------------------------------- */
/*		PROCEDURE FOR A >= 8 AND B >= 8 */
/* ----------------------------------------------------------------------- */
	double h, x0, y0, lambda;
	if (a <= b) {
	    h = a / b;
	    x0 = h / (h + 1.);
	    y0 = 1. / (h + 1.);
	    lambda = a - (a + b) * x;
	} else {
	    h = b / a;
	    x0 = 1. / (h + 1.);
	    y0 = h / (h + 1.);
	    lambda = (a + b) * y - b;
	}

	e = -lambda / a;
	if (fabs(e) > .6)
	    u = e - log(x / x0);
	else
	    u = rlog1(e);

	e = lambda / b;
	if (fabs(e) <= .6)
	    v = rlog1(e);
	else
	    v = e - log(y / y0);

	z = log_p ? -(a * u + b * v) : exp(-(a * u + b * v));

	return(log_p
	       ? -M_LN_SQRT_2PI + .5*log(b * x0) + z - bcorr(a,b)
	       : const__ * sqrt(b * x0) * z * exp(-bcorr(a, b)));
    }
} /* brcomp */

// called only once from  bup(),  as   r = brcmp1(mu, a, b, x, y, FALSE) / a;
//                        -----
static double brcmp1(int mu, double a, double b, double x, double y, int give_log)
{
/* -----------------------------------------------------------------------
 *          Evaluation of    exp(mu) * x^a * y^b / beta(a,b)
 * ----------------------------------------------------------------------- */

    static double const__ = .398942280401433; /* == 1/sqrt(2*pi); */
    /* R has  M_1_SQRT_2PI */

    /* Local variables */
    double c, t, u, v, z, a0, b0, apb;

    a0 = min(a,b);
    if (a0 < 8.) {
	double lnx, lny;
	if (x <= .375) {
	    lnx = log(x);
	    lny = alnrel(-x);
	} else if (y > .375) {
	    // L11:
	    lnx = log(x);
	    lny = log(y);
	} else {
	    lnx = alnrel(-y);
	    lny = log(y);
	}

	// L20:
	z = a * lnx + b * lny;
	if (a0 >= 1.) {
	    z -= betaln(a, b);
	    return esum(mu, z, give_log);
	}
	// else :
	/* ----------------------------------------------------------------------- */
	/*              PROCEDURE FOR A < 1 OR B < 1 */
	/* ----------------------------------------------------------------------- */
	// L30:
	b0 = max(a,b);
	if (b0 >= 8.) {
	/* L80:                  ALGORITHM FOR b0 >= 8 */
	    u = gamln1(a0) + algdiv(a0, b0);
	    R_ifDEBUG_printf(" brcmp1(mu,a,b,*): a0 < 1, b0 >= 8;  z=%.15g\n", z);
	    return give_log
		? log(a0) + esum(mu, z - u, TRUE)
		:     a0  * esum(mu, z - u, FALSE);

	} else if (b0 <= 1.) {
	    //                   a0 < 1, b0 <= 1
	    double ans = esum(mu, z, give_log);
	    if (ans == (give_log ? ML_NEGINF : 0.))
		return ans;

	    apb = a + b;
	    if (apb > 1.) {
		// L40:
		u = a + b - 1.;
		z = (gam1(u) + 1.) / apb;
	    } else {
		z = gam1(apb) + 1.;
	    }
	    // L50:
	    c = give_log
		? log1p(gam1(a)) + log1p(gam1(b)) - log(z)
		: (gam1(a) + 1.) * (gam1(b) + 1.) / z;
	    R_ifDEBUG_printf(" brcmp1(mu,a,b,*): a0 < 1, b0 <= 1;  c=%.15g\n", c);
	    return give_log
		? ans + log(a0) + c - log1p(a0 / b0)
		: ans * (a0 * c) / (a0 / b0 + 1.);
	}
	// else:               algorithm for	a0 < 1 < b0 < 8
	// L60:
	u = gamln1(a0);
	int n = (int)(b0 - 1.);
	if (n >= 1) {
	    c = 1.;
	    for (int i = 1; i <= n; ++i) {
		b0 += -1.;
		c *= b0 / (a0 + b0);
		/* L61: */
	    }
	    u += log(c); // TODO?: log(c) = log( prod(...) ) =  sum( log(...) )
	}
	// L70:
	z -= u;
	b0 += -1.;
	apb = a0 + b0;
	if (apb > 1.) {
	    // L71:
	    t = (gam1(apb - 1.) + 1.) / apb;
	} else {
	    t = gam1(apb) + 1.;
	}
	R_ifDEBUG_printf(" brcmp1(mu,a,b,*): a0 < 1 < b0 < 8;  t=%.15g\n", t);
	// L72:
	return give_log
	    ? log(a0)+ esum(mu, z, TRUE) + log1p(gam1(b0)) - log(t) // TODO? log(t) = log1p(..)
	    :     a0 * esum(mu, z, FALSE) * (gam1(b0) + 1.) / t;

    } else {

/* ----------------------------------------------------------------------- */
/*              PROCEDURE FOR A >= 8 AND B >= 8 */
/* ----------------------------------------------------------------------- */
	// L100:
	double h, x0, y0, lambda;
	if (a > b) {
	    // L101:
	    h = b / a;
	    x0 = 1. / (h + 1.);// => lx0 := log(x0) = 0 - log1p(h)
	    y0 = h / (h + 1.);
	    lambda = (a + b) * y - b;
	} else {
	    h = a / b;
	    x0 = h / (h + 1.);  // => lx0 := log(x0) = - log1p(1/h)
	    y0 = 1. / (h + 1.);
	    lambda = a - (a + b) * x;
	}
	double lx0 = -log1p(b/a); // in both cases

	R_ifDEBUG_printf(" brcmp1(mu,a,b,*): a,b >= 8;	x0=%.15g, lx0=log(x0)=%.15g\n",
			 x0, lx0);
	// L110:
	double e = -lambda / a;
	if (fabs(e) > 0.6) {
	    // L111:
	    u = e - log(x / x0);
	} else {
	    u = rlog1(e);
	}

	// L120:
	e = lambda / b;
	if (fabs(e) > 0.6) {
	    // L121:
	    v = e - log(y / y0);
	} else {
	    v = rlog1(e);
	}

	// L130:
	z = esum(mu, -(a * u + b * v), give_log);
	return give_log
	    ? log(const__)+ (log(b) + lx0)/2. + z      - bcorr(a, b)
	    :     const__ * sqrt(b * x0)      * z * exp(-bcorr(a, b));
    }

} /* brcmp1 */

static void bgrat(double a, double b, double x, double y, double *w,
		  double eps, int *ierr, Rboolean log_w)
{
/* -----------------------------------------------------------------------
*     Asymptotic Expansion for I_x(a,b)  when a is larger than b.
*     Compute   w := w + I_x(a,b)
*     It is assumed a >= 15 and b <= 1.
*     eps is the tolerance used.
*     ierr is a variable that reports the status of the results.
*
* if(log_w),  *w  itself must be in log-space;
*     compute   w := w + I_x(a,b)  but return *w = log(w):
*          *w := log(exp(*w) + I_x(a,b)) = logspace_add(*w, log( I_x(a,b) ))
* ----------------------------------------------------------------------- */

#define n_terms_bgrat 30
    double c[n_terms_bgrat], d[n_terms_bgrat];
    double bm1 = b - 0.5 - 0.5,
	nu = a + bm1 * 0.5, /* nu = a + (b-1)/2 =: T, in (9.1) of
			     * Didonato & Morris(1992), p.362 */
	lnx = (y > 0.375) ? log(x) : alnrel(-y),
	z = -nu * lnx; // z =: u in (9.1) of D.&M.(1992)

    if (b * z == 0.) { // should not happen, but does, e.g.,
	// for  pbeta(1e-320, 1e-5, 0.5)  i.e., _subnormal_ x,
	// Warning ... bgrat(a=20.5, b=1e-05, x=1, y=9.99989e-321): ..
	MATHLIB_WARNING5(
	    "bgrat(a=%g, b=%g, x=%g, y=%g): z=%g, b*z == 0 underflow, hence inaccurate pbeta()",
	    a,b,x,y, z);
	/* L_Error:    THE EXPANSION CANNOT BE COMPUTED */
	 *ierr = 1; return;
    }

/*                 COMPUTATION OF THE EXPANSION */
    double
	/* r1 = b * (gam1(b) + 1.) * exp(b * log(z)),// = b/gamma(b+1) z^b = z^b / gamma(b)
	 * set r := exp(-z) * z^b / gamma(b) ;
	 *          gam1(b) = 1/gamma(b+1) - 1 , b in [-1/2, 3/2] */
	// exp(a*lnx) underflows for large (a * lnx); e.g. large a ==> using log_r := log(r):
	// r = r1 * exp(a * lnx) * exp(bm1 * 0.5 * lnx);
	// log(r)=log(b) + log1p(gam1(b)) + b * log(z) + (a * lnx) + (bm1 * 0.5 * lnx),
	log_r = log(b) + log1p(gam1(b)) + b * log(z) + nu * lnx,
	// FIXME work with  log_u = log(u)  also when log_p=FALSE  (??)
	// u is 'factored out' from the expansion {and multiplied back, at the end}:
	log_u = log_r - (algdiv(b, a) + b * log(nu)),// algdiv(b,a) = log(gamma(a)/gamma(a+b))
	/* u = (log_p) ? log_r - u : exp(log_r-u); // =: M  in (9.2) of {reference above} */
	/* u = algdiv(b, a) + b * log(nu);// algdiv(b,a) = log(gamma(a)/gamma(a+b)) */
	// u = (log_p) ? log_u : exp(log_u); // =: M  in (9.2) of {reference above}
	u = exp(log_u);

    if (log_u == ML_NEGINF) {
	R_ifDEBUG_printf(" bgrat(*): underflow log_u = -Inf  = log_r -u', log_r = %g ",
			 log_r);
	/* L_Error:    THE EXPANSION CANNOT BE COMPUTED */ *ierr = 2; return;
    }

    Rboolean u_0 = (u == 0.); // underflow --> do work with log(u) == log_u !
    double l = // := *w/u .. but with care: such that it also works when u underflows to 0:
	log_w
	? ((*w == ML_NEGINF) ? 0. : exp(  *w    - log_u))
	: ((*w == 0.)        ? 0. : exp(log(*w) - log_u));

    R_ifDEBUG_printf(" bgrat(a=%g, b=%g, x=%g, *)\n -> u=%g, l='w/u'=%g, ",
		     a,b,x, u, l);
    double
	q_r = grat_r(b, z, log_r, eps), // = q/r of former grat1(b,z, r, &p, &q)
	v = 0.25 / (nu * nu),
	t2 = lnx * 0.25 * lnx,
	j = q_r,
	sum = j,
	t = 1., cn = 1., n2 = 0.;
    for (int n = 1; n <= n_terms_bgrat; ++n) {
	double bp2n = b + n2;
	j = (bp2n * (bp2n + 1.) * j + (z + bp2n + 1.) * t) * v;
	n2 += 2.;
	t *= t2;
	cn /= n2 * (n2 + 1.);
	int nm1 = n - 1;
	c[nm1] = cn;
	double s = 0.;
	if (n > 1) {
	    double coef = b - n;
	    for (int i = 1; i <= nm1; ++i) {
		s += coef * c[i - 1] * d[nm1 - i];
		coef += b;
	    }
	}
	d[nm1] = bm1 * cn + s / n;
	double dj = d[nm1] * j;
	sum += dj;
	if (sum <= 0.) {
	    R_ifDEBUG_printf(" bgrat(*): sum_n(..) <= 0; should not happen (n=%d)\n", n);
	    /* L_Error:    THE EXPANSION CANNOT BE COMPUTED */ *ierr = 3; return;
	}
	if (fabs(dj) <= eps * (sum + l)) {
	    *ierr = 0;
	    break;
	} else if(n == n_terms_bgrat) { // never? ; please notify R-core if seen:
	    *ierr = 4;
	    MATHLIB_WARNING5(
	"bgrat(a=%g, b=%g, x=%g) *no* convergence: NOTIFY R-core!\n dj=%g, rel.err=%g\n",
		a,b,x, dj, fabs(dj) /(sum + l));
	}
    } // for(n .. n_terms..)

/*                    ADD THE RESULTS TO W */

    if(log_w) // *w is in log space already:
	*w = logspace_add(*w, log_u + log(sum));
    else
	*w += (u_0 ? exp(log_u + log(sum)) : u * sum);
    return;
} /* bgrat */


// called only from bgrat() , as   q_r = grat_r(b, z, log_r, eps)  :
static double grat_r(double a, double x, double log_r, double eps)
{
/* -----------------------------------------------------------------------
 *        Scaled complement of incomplete gamma ratio function
 *                   grat_r(a,x,r) :=  Q(a,x) / r
 * where
 *               Q(a,x) = pgamma(x,a, lower.tail=FALSE)
 *     and            r = e^(-x)* x^a / Gamma(a) ==  exp(log_r)
 *
 *     It is assumed that a <= 1.  eps is the tolerance to be used.
 * ----------------------------------------------------------------------- */

    if (a * x == 0.) { /* L130: */
	if (x <= a) {
	    /* L100: */ return exp(-log_r);
	} else {
	    /* L110:*/  return 0.;
	}
    }
    else if (a == 0.5) { // e.g. when called from pt()
	/* L120: */
	if (x < 0.25) {
	    double p = erf__(sqrt(x));
	    R_ifDEBUG_printf(" grat_r(a=%g, x=%g ..)): a=1/2 --> p=erf__(.)= %g\n",
			     a, x, p);
	    return (0.5 - p + 0.5)*exp(-log_r);

        } else { // 2013-02-27: improvement for "large" x: direct computation of q/r:
	    double sx = sqrt(x),
		q_r = erfc1(1, sx)/sx * M_SQRT_PI;
	    R_ifDEBUG_printf(" grat_r(a=%g, x=%g ..)): a=1/2 --> q_r=erfc1(..)/r= %g\n",
			     a,x, q_r);
	    return q_r;
	}

    } else if (x < 1.1) { /* L10:  Taylor series for  P(a,x)/x^a */

	double an = 3.,
	    c = x,
	    sum = x / (a + 3.),
	    tol = eps * 0.1 / (a + 1.), t;
	do {
	    an += 1.;
	    c *= -(x / an);
	    t = c / (a + an);
	    sum += t;
	} while (fabs(t) > tol);

	R_ifDEBUG_printf(" grat_r(a=%g, x=%g, log_r=%g): sum=%g; Taylor w/ %.0f terms",
			 a,x,log_r, sum, an-3.);
	double j = a * x * ((sum/6. - 0.5/(a + 2.)) * x + 1./(a + 1.)),
	    z = a * log(x),
	    h = gam1(a),
	    g = h + 1.;

	if ((x >= 0.25 && (a < x / 2.59)) || (z > -0.13394)) {
	    // L40:
	    double l = rexpm1(z),
		q = ((l + 0.5 + 0.5) * j - l) * g - h;
	    if (q <= 0.) {
		R_ifDEBUG_printf(" => q_r= 0.\n");
		/* L110:*/ return 0.;
	    } else {
		R_ifDEBUG_printf(" => q_r=%.15g\n", q * exp(-log_r));
		return q * exp(-log_r);
	    }

	} else {
	    double p = exp(z) * g * (0.5 - j + 0.5);
	    R_ifDEBUG_printf(" => q_r=%.15g\n", (0.5 - p + 0.5) * exp(-log_r));
	    return /* q/r = */ (0.5 - p + 0.5) * exp(-log_r);
	}

    } else {
	/* L50: ----  (x >= 1.1)  ---- Continued Fraction Expansion */

	double a2n_1 = 1.,
	    a2n = 1.,
	    b2n_1 = x,
	    b2n = x + (1. - a),
	    c = 1., am0, an0;

	do {
	    a2n_1 = x * a2n + c * a2n_1;
	    b2n_1 = x * b2n + c * b2n_1;
	    am0 = a2n_1 / b2n_1;
	    c += 1.;
	    double c_a = c - a;
	    a2n = a2n_1 + c_a * a2n;
	    b2n = b2n_1 + c_a * b2n;
	    an0 = a2n / b2n;
	} while (fabs(an0 - am0) >= eps * an0);

	R_ifDEBUG_printf(" grat_r(a=%g, x=%g, log_r=%g): Cont.frac. %.0f terms => q_r=%.15g\n",
			 a,x, log_r, c-1., an0);
	return /* q/r = (r * an0)/r = */ an0;
    }
} /* grat_r */



static double basym(double a, double b, double lambda, double eps, int log_p)
{
/* ----------------------------------------------------------------------- */
/*     ASYMPTOTIC EXPANSION FOR I_x(A,B) FOR LARGE A AND B. */
/*     LAMBDA = (A + B)*Y - B  AND EPS IS THE TOLERANCE USED. */
/*     IT IS ASSUMED THAT LAMBDA IS NONNEGATIVE AND THAT */
/*     A AND B ARE GREATER THAN OR EQUAL TO 15. */
/* ----------------------------------------------------------------------- */


/* ------------------------ */
/*     ****** NUM IS THE MAXIMUM VALUE THAT N CAN TAKE IN THE DO LOOP */
/*            ENDING AT STATEMENT 50. IT IS REQUIRED THAT NUM BE EVEN. */
#define num_IT 20
/*            THE ARRAYS A0, B0, C, D HAVE DIMENSION NUM + 1. */

    static double const e0 = 1.12837916709551;/* e0 == 2/sqrt(pi) */
    static double const e1 = .353553390593274;/* e1 == 2^(-3/2)   */
    static double const ln_e0 = 0.120782237635245; /* == ln(e0) */

    double a0[num_IT + 1], b0[num_IT + 1], c[num_IT + 1], d[num_IT + 1];

    double f = a * rlog1(-lambda/a) + b * rlog1(lambda/b), t;
    if(log_p)
	t = -f;
    else {
	t = exp(-f);
	if (t == 0.) {
	    return 0; /* once underflow, always underflow .. */
	}
    }
    double z0 = sqrt(f),
	z = z0 / e1 * 0.5,
	z2 = f + f,
	h, r0, r1, w0;

    if (a < b) {
	h = a / b;
	r0 = 1. / (h + 1.);
	r1 = (b - a) / b;
	w0 = 1. / sqrt(a * (h + 1.));
    } else {
	h = b / a;
	r0 = 1. / (h + 1.);
	r1 = (b - a) / a;
	w0 = 1. / sqrt(b * (h + 1.));
    }

    a0[0] = r1 * .66666666666666663;
    c[0] = a0[0] * -0.5;
    d[0] = -c[0];
    double j0 = 0.5 / e0 * erfc1(1, z0),
	j1 = e1,
	sum = j0 + d[0] * w0 * j1;

    double s = 1.,
	h2 = h * h,
	hn = 1.,
	w = w0,
	znm1 = z,
	zn = z2;
    for (int n = 2; n <= num_IT; n += 2) {
	hn *= h2;
	a0[n - 1] = r0 * 2. * (h * hn + 1.) / (n + 2.);
	int np1 = n + 1;
	s += hn;
	a0[np1 - 1] = r1 * 2. * s / (n + 3.);

	for (int i = n; i <= np1; ++i) {
	    double r = (i + 1.) * -0.5;
	    b0[0] = r * a0[0];
	    for (int m = 2; m <= i; ++m) {
		double bsum = 0.;
		for (int j = 1; j <= m-1; ++j) {
		    int mmj = m - j;
		    bsum += (j * r - mmj) * a0[j - 1] * b0[mmj - 1];
		}
		b0[m - 1] = r * a0[m - 1] + bsum / m;
	    }
	    c[i - 1] = b0[i - 1] / (i + 1.);

	    double dsum = 0.;
	    for (int j = 1; j <= i-1; ++j) {
		dsum += d[i - j - 1] * c[j - 1];
	    }
	    d[i - 1] = -(dsum + c[i - 1]);
	}

	j0 = e1 * znm1 + (n - 1.) * j0;
	j1 = e1 * zn + n * j1;
	znm1 = z2 * znm1;
	zn = z2 * zn;
	w *= w0;
	double t0 = d[n - 1] * w * j0;
	w *= w0;
	double t1 = d[np1 - 1] * w * j1;
	sum += t0 + t1;
	if (fabs(t0) + fabs(t1) <= eps * sum) {
	    break;
	}
    }

    if(log_p)
	return ln_e0 + t - bcorr(a, b) + log(sum);
    else {
	double u = exp(-bcorr(a, b));
	return e0 * t * u * sum;
    }

} /* basym_ */


static double exparg(int l)
{
/* --------------------------------------------------------------------
 *     If l = 0 then  exparg(l) = The largest positive W for which
 *     exp(W) can be computed. With 0.99999 fuzz  ==> exparg(0) =   709.7756  nowadays

 *     if l = 1 (nonzero) then  exparg(l) = the largest negative W for
 *     which the computed value of exp(W) is nonzero.
 *     With 0.99999 fuzz			  ==> exparg(1) =  -709.0825  nowadays

 *     Note... only an approximate value for exparg(L) is needed.
 * -------------------------------------------------------------------- */

    static double const lnb = .69314718055995;
    int m = (l == 0) ? Rf_i1mach(16) : Rf_i1mach(15) - 1;

    return m * lnb * .99999;
} /* exparg */

static double esum(int mu, double x, int give_log)
{
/* ----------------------------------------------------------------------- */
/*                    EVALUATION OF EXP(MU + X) */
/* ----------------------------------------------------------------------- */

    if(give_log)
	return x + (double) mu;

    // else :
    double w;
    if (x > 0.) { /* L10: */
	if (mu > 0)  return exp((double) mu) * exp(x);
	w = mu + x;
	if (w < 0.) return exp((double) mu) * exp(x);
    }
    else { /* x <= 0 */
	if (mu < 0)  return exp((double) mu) * exp(x);
	w = mu + x;
	if (w > 0.) return exp((double) mu) * exp(x);
    }
    return exp(w);

} /* esum */

double rexpm1(double x)
{
/* ----------------------------------------------------------------------- */
/*            EVALUATION OF THE FUNCTION EXP(X) - 1 */
/* ----------------------------------------------------------------------- */

    static double p1 = 9.14041914819518e-10;
    static double p2 = .0238082361044469;
    static double q1 = -.499999999085958;
    static double q2 = .107141568980644;
    static double q3 = -.0119041179760821;
    static double q4 = 5.95130811860248e-4;

    if (fabs(x) <= 0.15) {
	return x * (((p2 * x + p1) * x + 1.) /
		    ((((q4 * x + q3) * x + q2) * x + q1) * x + 1.));
    }
    else { /* |x| > 0.15 : */
	double w = exp(x);
	if (x > 0.)
	    return w * (0.5 - 1. / w + 0.5);
	else
	    return w - 0.5 - 0.5;
    }

} /* rexpm1 */

static double alnrel(double a)
{
/* -----------------------------------------------------------------------
 *            Evaluation of the function ln(1 + a)
 * ----------------------------------------------------------------------- */

    if (fabs(a) > 0.375)
	return log(1. + a);
    // else : |a| <= 0.375
    static double
	p1 = -1.29418923021993,
	p2 = .405303492862024,
	p3 = -.0178874546012214,
	q1 = -1.62752256355323,
	q2 = .747811014037616,
	q3 = -.0845104217945565;
    double
	t = a / (a + 2.),
	t2 = t * t,
	w = (((p3 * t2 + p2) * t2 + p1) * t2 + 1.) /
	(((q3 * t2 + q2) * t2 + q1) * t2 + 1.);
    return t * 2. * w;

} /* alnrel */

static double rlog1(double x)
{
/* -----------------------------------------------------------------------
 *             Evaluation of the function  x - ln(1 + x)
 * ----------------------------------------------------------------------- */

    static double a = .0566749439387324;
    static double b = .0456512608815524;
    static double p0 = .333333333333333;
    static double p1 = -.224696413112536;
    static double p2 = .00620886815375787;
    static double q1 = -1.27408923933623;
    static double q2 = .354508718369557;

    double h, r, t, w, w1;
    if (x < -0.39 || x > 0.57) { /* direct evaluation */
	w = x + 0.5 + 0.5;
	return x - log(w);
    }
    /* else */
    if (x < -0.18) { /* L10: */
	h = x + .3;
	h /= .7;
	w1 = a - h * .3;
    }
    else if (x > 0.18) { /* L20: */
	h = x * .75 - .25;
	w1 = b + h / 3.;
    }
    else { /*		Argument Reduction */
	h = x;
	w1 = 0.;
    }

/* L30:              	Series Expansion */

    r = h / (h + 2.);
    t = r * r;
    w = ((p2 * t + p1) * t + p0) / ((q2 * t + q1) * t + 1.);
    return t * 2. * (1. / (1. - r) - r * w) + w1;

} /* rlog1 */

static double erf__(double x)
{
/* -----------------------------------------------------------------------
 *             EVALUATION OF THE REAL ERROR FUNCTION
 * ----------------------------------------------------------------------- */

    /* Initialized data */

    static double c = .564189583547756;
    static double a[5] = { 7.7105849500132e-5,-.00133733772997339,
	    .0323076579225834,.0479137145607681,.128379167095513 };
    static double b[3] = { .00301048631703895,.0538971687740286,
	    .375795757275549 };
    static double p[8] = { -1.36864857382717e-7,.564195517478974,
	    7.21175825088309,43.1622272220567,152.98928504694,
	    339.320816734344,451.918953711873,300.459261020162 };
    static double q[8] = { 1.,12.7827273196294,77.0001529352295,
	    277.585444743988,638.980264465631,931.35409485061,
	    790.950925327898,300.459260956983 };
    static double r[5] = { 2.10144126479064,26.2370141675169,
	    21.3688200555087,4.6580782871847,.282094791773523 };
    static double s[4] = { 94.153775055546,187.11481179959,
	    99.0191814623914,18.0124575948747 };

    /* Local variables */
    double t, x2, ax, bot, top;

    ax = fabs(x);
    if (ax <= 0.5) {
	t = x * x;
	top = (((a[0] * t + a[1]) * t + a[2]) * t + a[3]) * t + a[4] + 1.;
	bot = ((b[0] * t + b[1]) * t + b[2]) * t + 1.;

	return x * (top / bot);
    }

    // else:  |x| > 0.5

    if (ax <= 4.) { //  |x| in (0.5, 4]
	top = ((((((p[0] * ax + p[1]) * ax + p[2]) * ax + p[3]) * ax + p[4]) * ax
		+ p[5]) * ax + p[6]) * ax + p[7];
	bot = ((((((q[0] * ax + q[1]) * ax + q[2]) * ax + q[3]) * ax + q[4]) * ax
		+ q[5]) * ax + q[6]) * ax + q[7];
	double R = 0.5 - exp(-x * x) * top / bot + 0.5;
	return (x < 0) ? -R : R;
    }

    // else:  |x| > 4

    if (ax >= 5.8) {
	return x > 0 ? 1 : -1;
    }

    // else:  4 < |x| < 5.8
    x2 = x * x;
    t = 1. / x2;
    top = (((r[0] * t + r[1]) * t + r[2]) * t + r[3]) * t + r[4];
    bot = (((s[0] * t + s[1]) * t + s[2]) * t + s[3]) * t + 1.;
    t = (c - top / (x2 * bot)) / ax;
    double R = 0.5 - exp(-x2) * t + 0.5;
    return (x < 0) ? -R : R;
} /* erf */

static double erfc1(int ind, double x)
{
/* ----------------------------------------------------------------------- */
/*         EVALUATION OF THE COMPLEMENTARY ERROR FUNCTION */

/*          ERFC1(IND,X) = ERFC(X)            IF IND = 0 */
/*          ERFC1(IND,X) = EXP(X*X)*ERFC(X)   OTHERWISE */
/* ----------------------------------------------------------------------- */

    /* Initialized data */

    static double c = .564189583547756;
    static double a[5] = { 7.7105849500132e-5,-.00133733772997339,
	    .0323076579225834,.0479137145607681,.128379167095513 };
    static double b[3] = { .00301048631703895,.0538971687740286,
	    .375795757275549 };
    static double p[8] = { -1.36864857382717e-7,.564195517478974,
	    7.21175825088309,43.1622272220567,152.98928504694,
	    339.320816734344,451.918953711873,300.459261020162 };
    static double q[8] = { 1.,12.7827273196294,77.0001529352295,
	    277.585444743988,638.980264465631,931.35409485061,
	    790.950925327898,300.459260956983 };
    static double r[5] = { 2.10144126479064,26.2370141675169,
	    21.3688200555087,4.6580782871847,.282094791773523 };
    static double s[4] = { 94.153775055546,187.11481179959,
	    99.0191814623914,18.0124575948747 };

    double ret_val;
    double e, t, w, bot, top;

    double ax = fabs(x);
    //				|X| <= 0.5 */
    if (ax <= 0.5) {
	double t = x * x,
	    top = (((a[0] * t + a[1]) * t + a[2]) * t + a[3]) * t + a[4] + 1.,
	    bot = ((b[0] * t + b[1]) * t + b[2]) * t + 1.;
	ret_val = 0.5 - x * (top / bot) + 0.5;
	if (ind != 0) {
	    ret_val = exp(t) * ret_val;
	}
	return ret_val;
    }
    // else (L10:):		0.5 < |X| <= 4
    if (ax <= 4.) {
	top = ((((((p[0] * ax + p[1]) * ax + p[2]) * ax + p[3]) * ax + p[4]) * ax
		+ p[5]) * ax + p[6]) * ax + p[7];
	bot = ((((((q[0] * ax + q[1]) * ax + q[2]) * ax + q[3]) * ax + q[4]) * ax
		+ q[5]) * ax + q[6]) * ax + q[7];
	ret_val = top / bot;

    } else { //			|X| > 4
	// L20:
	if (x <= -5.6) {
	    // L50:            	LIMIT VALUE FOR "LARGE" NEGATIVE X
	    ret_val = 2.;
	    if (ind != 0) {
		ret_val = exp(x * x) * 2.;
	    }
	    return ret_val;
	}
	if (ind == 0 && (x > 100. || x * x > -exparg(1))) {
	    // LIMIT VALUE FOR LARGE POSITIVE X   WHEN IND = 0
	    // L60:
	    return 0.;
	}

	// L30:
	t = 1. / (x * x);
	top = (((r[0] * t + r[1]) * t + r[2]) * t + r[3]) * t + r[4];
	bot = (((s[0] * t + s[1]) * t + s[2]) * t + s[3]) * t + 1.;
	ret_val = (c - t * top / bot) / ax;
    }

    // L40:                 FINAL ASSEMBLY
    if (ind != 0) {
	if (x < 0.)
	    ret_val = exp(x * x) * 2. - ret_val;
    } else {
	// L41:  ind == 0 :
	w = x * x;
	t = w;
	e = w - t;
	ret_val = (0.5 - e + 0.5) * exp(-t) * ret_val;
	if (x < 0.)
	    ret_val = 2. - ret_val;
    }
    return ret_val;

} /* erfc1 */

static double gam1(double a)
{
/*     ------------------------------------------------------------------ */
/*     COMPUTATION OF 1/GAMMA(A+1) - 1  FOR -0.5 <= A <= 1.5 */
/*     ------------------------------------------------------------------ */

    double d, t, w, bot, top;

    t = a;
    d = a - 0.5;
    // t := if(a > 1/2)  a-1  else  a
    if (d > 0.)
	t = d - 0.5;
    if (t < 0.) { /* L30: */
	static double
	    r[9] = { -.422784335098468,-.771330383816272,
		     -.244757765222226,.118378989872749,9.30357293360349e-4,
		     -.0118290993445146,.00223047661158249,2.66505979058923e-4,
		     -1.32674909766242e-4 },
	    s1 = .273076135303957,
	    s2 = .0559398236957378;

	top = (((((((r[8] * t + r[7]) * t + r[6]) * t + r[5]) * t + r[4]
		     ) * t + r[3]) * t + r[2]) * t + r[1]) * t + r[0];
	bot = (s2 * t + s1) * t + 1.;
	w = top / bot;
	R_ifDEBUG_printf("  gam1(a = %.15g): t < 0: w=%.15g\n", a, w);
	if (d > 0.)
	    return t * w / a;
	else
	    return a * (w + 0.5 + 0.5);

    } else if (t == 0) { // L10: a in {0, 1}
	return 0.;

    } else { /* t > 0;  L20: */
	static double
	    p[7] = { .577215664901533,-.409078193005776,
		     -.230975380857675,.0597275330452234,.0076696818164949,
		     -.00514889771323592,5.89597428611429e-4 },
	    q[5] = { 1.,.427569613095214,.158451672430138,
		     .0261132021441447,.00423244297896961 };

	top = (((((p[6] * t + p[5]) * t + p[4]) * t + p[3]) * t + p[2]
		   ) * t + p[1]) * t + p[0];
	bot = (((q[4] * t + q[3]) * t + q[2]) * t + q[1]) * t + 1.;
	w = top / bot;
	R_ifDEBUG_printf("  gam1(a = %.15g): t > 0: (is a < 1.5 ?)  w=%.15g\n",
			 a, w);
	if (d > 0.) /* L21: */
	    return t / a * (w - 0.5 - 0.5);
	else
	    return a * w;
    }
} /* gam1 */

static double gamln1(double a)
{
/* ----------------------------------------------------------------------- */
/*     EVALUATION OF LN(GAMMA(1 + A)) FOR -0.2 <= A <= 1.25 */
/* ----------------------------------------------------------------------- */

    double w;
    if (a < 0.6) {
	static double p0 = .577215664901533;
	static double p1 = .844203922187225;
	static double p2 = -.168860593646662;
	static double p3 = -.780427615533591;
	static double p4 = -.402055799310489;
	static double p5 = -.0673562214325671;
	static double p6 = -.00271935708322958;
	static double q1 = 2.88743195473681;
	static double q2 = 3.12755088914843;
	static double q3 = 1.56875193295039;
	static double q4 = .361951990101499;
	static double q5 = .0325038868253937;
	static double q6 = 6.67465618796164e-4;
	w = ((((((p6 * a + p5)* a + p4)* a + p3)* a + p2)* a + p1)* a + p0) /
	    ((((((q6 * a + q5)* a + q4)* a + q3)* a + q2)* a + q1)* a + 1.);
	return -(a) * w;
    }
    else { /* 0.6 <= a <= 1.25 */
	static double r0 = .422784335098467;
	static double r1 = .848044614534529;
	static double r2 = .565221050691933;
	static double r3 = .156513060486551;
	static double r4 = .017050248402265;
	static double r5 = 4.97958207639485e-4;
	static double s1 = 1.24313399877507;
	static double s2 = .548042109832463;
	static double s3 = .10155218743983;
	static double s4 = .00713309612391;
	static double s5 = 1.16165475989616e-4;
	double x = a - 0.5 - 0.5;
	w = (((((r5 * x + r4) * x + r3) * x + r2) * x + r1) * x + r0) /
	    (((((s5 * x + s4) * x + s3) * x + s2) * x + s1) * x + 1.);
	return x * w;
    }
} /* gamln1 */

static double psi(double x)
{
/* ---------------------------------------------------------------------

 *                 Evaluation of the Digamma function psi(x)

 *                           -----------

 *     Psi(xx) is assigned the value 0 when the digamma function cannot
 *     be computed.

 *     The main computation involves evaluation of rational Chebyshev
 *     approximations published in Math. Comp. 27, 123-127(1973) by
 *     Cody, Strecok and Thacher. */

/* --------------------------------------------------------------------- */
/*     Psi was written at Argonne National Laboratory for the FUNPACK */
/*     package of special function subroutines. Psi was modified by */
/*     A.H. Morris (NSWC). */
/* --------------------------------------------------------------------- */

    static double piov4 = .785398163397448; /* == pi / 4 */
/*     dx0 = zero of psi() to extended precision : */
    static double dx0 = 1.461632144968362341262659542325721325;

/* --------------------------------------------------------------------- */
/*     COEFFICIENTS FOR RATIONAL APPROXIMATION OF */
/*     PSI(X) / (X - X0),  0.5 <= X <= 3. */
    static double p1[7] = { .0089538502298197,4.77762828042627,
	    142.441585084029,1186.45200713425,3633.51846806499,
	    4138.10161269013,1305.60269827897 };
    static double q1[6] = { 44.8452573429826,520.752771467162,
	    2210.0079924783,3641.27349079381,1908.310765963,
	    6.91091682714533e-6 };
/* --------------------------------------------------------------------- */


/* --------------------------------------------------------------------- */
/*     COEFFICIENTS FOR RATIONAL APPROXIMATION OF */
/*     PSI(X) - LN(X) + 1 / (2*X),  X > 3. */

    static double p2[4] = { -2.12940445131011,-7.01677227766759,
	    -4.48616543918019,-.648157123766197 };
    static double q2[4] = { 32.2703493791143,89.2920700481861,
	    54.6117738103215,7.77788548522962 };
/* --------------------------------------------------------------------- */

    int i, m, n, nq;
    double d2;
    double w, z;
    double den, aug, sgn, xmx0, xmax1, upper, xsmall;

/* --------------------------------------------------------------------- */


/*     MACHINE DEPENDENT CONSTANTS ... */

/* --------------------------------------------------------------------- */
/*	  XMAX1	 = THE SMALLEST POSITIVE FLOATING POINT CONSTANT
		   WITH ENTIRELY INT REPRESENTATION.  ALSO USED
		   AS NEGATIVE OF LOWER BOUND ON ACCEPTABLE NEGATIVE
		   ARGUMENTS AND AS THE POSITIVE ARGUMENT BEYOND WHICH
		   PSI MAY BE REPRESENTED AS LOG(X).
 * Originally:  xmax1 = amin1(ipmpar(3), 1./spmpar(1))  */
    xmax1 = (double) INT_MAX;
    d2 = 0.5 / Rf_d1mach(3); /*= 0.5 / (0.5 * DBL_EPS) = 1/DBL_EPSILON = 2^52 */
    if(xmax1 > d2) xmax1 = d2;

/* --------------------------------------------------------------------- */
/*        XSMALL = ABSOLUTE ARGUMENT BELOW WHICH PI*COTAN(PI*X) */
/*                 MAY BE REPRESENTED BY 1/X. */
    xsmall = 1e-9;
/* --------------------------------------------------------------------- */
    aug = 0.;
    if (x < 0.5) {
/* --------------------------------------------------------------------- */
/*     X < 0.5,  USE REFLECTION FORMULA */
/*     PSI(1-X) = PSI(X) + PI * COTAN(PI*X) */
/* --------------------------------------------------------------------- */
	if (fabs(x) <= xsmall) {

	    if (x == 0.) {
		goto L_err;
	    }
/* --------------------------------------------------------------------- */
/*     0 < |X| <= XSMALL.  USE 1/X AS A SUBSTITUTE */
/*     FOR  PI*COTAN(PI*X) */
/* --------------------------------------------------------------------- */
	    aug = -1. / x;
	} else { /* |x| > xsmall */
/* --------------------------------------------------------------------- */
/*     REDUCTION OF ARGUMENT FOR COTAN */
/* --------------------------------------------------------------------- */
	    /* L100: */
	    w = -x;
	    sgn = piov4;
	    if (w <= 0.) {
		w = -w;
		sgn = -sgn;
	    }
/* --------------------------------------------------------------------- */
/*     MAKE AN ERROR EXIT IF |X| >= XMAX1 */
/* --------------------------------------------------------------------- */
	    if (w >= xmax1) {
		goto L_err;
	    }
	    nq = (int) w;
	    w -= (double) nq;
	    nq = (int) (w * 4.);
	    w = (w - (double) nq * 0.25) * 4.;
/* --------------------------------------------------------------------- */
/*     W IS NOW RELATED TO THE FRACTIONAL PART OF  4. * X. */
/*     ADJUST ARGUMENT TO CORRESPOND TO VALUES IN FIRST */
/*     QUADRANT AND DETERMINE SIGN */
/* --------------------------------------------------------------------- */
	    n = nq / 2;
	    if (n + n != nq) {
		w = 1. - w;
	    }
	    z = piov4 * w;
	    m = n / 2;
	    if (m + m != n) {
		sgn = -sgn;
	    }
/* --------------------------------------------------------------------- */
/*     DETERMINE FINAL VALUE FOR  -PI*COTAN(PI*X) */
/* --------------------------------------------------------------------- */
	    n = (nq + 1) / 2;
	    m = n / 2;
	    m += m;
	    if (m == n) {
/* --------------------------------------------------------------------- */
/*     CHECK FOR SINGULARITY */
/* --------------------------------------------------------------------- */
		if (z == 0.) {
		    goto L_err;
		}
/* --------------------------------------------------------------------- */
/*     USE COS/SIN AS A SUBSTITUTE FOR COTAN, AND */
/*     SIN/COS AS A SUBSTITUTE FOR TAN */
/* --------------------------------------------------------------------- */
		aug = sgn * (cos(z) / sin(z) * 4.);

	    } else { /* L140: */
		aug = sgn * (sin(z) / cos(z) * 4.);
	    }
	}

	x = 1. - x;

    }
    /* L200: */
    if (x <= 3.) {
/* --------------------------------------------------------------------- */
/*     0.5 <= X <= 3. */
/* --------------------------------------------------------------------- */
	den = x;
	upper = p1[0] * x;

	for (i = 1; i <= 5; ++i) {
	    den = (den + q1[i - 1]) * x;
	    upper = (upper + p1[i]) * x;
	}

	den = (upper + p1[6]) / (den + q1[5]);
	xmx0 = x - dx0;
	return den * xmx0 + aug;
    }

/* --------------------------------------------------------------------- */
/*     IF X >= XMAX1, PSI = LN(X) */
/* --------------------------------------------------------------------- */
    if (x < xmax1) {
/* --------------------------------------------------------------------- */
/*     3. < X < XMAX1 */
/* --------------------------------------------------------------------- */
	w = 1. / (x * x);
	den = w;
	upper = p2[0] * w;

	for (i = 1; i <= 3; ++i) {
	    den = (den + q2[i - 1]) * w;
	    upper = (upper + p2[i]) * w;
	}

	aug = upper / (den + q2[3]) - 0.5 / x + aug;
    }
    return aug + log(x);

/* --------------------------------------------------------------------- */
/*     ERROR RETURN */
/* --------------------------------------------------------------------- */
L_err:
    return 0.;
} /* psi */

static double betaln(double a0, double b0)
{
/* -----------------------------------------------------------------------
 *     Evaluation of the logarithm of the beta function  ln(beta(a0,b0))
 * ----------------------------------------------------------------------- */

    static double e = .918938533204673;/* e == 0.5*LN(2*PI) */

    double
	a = min(a0 ,b0),
	b = max(a0, b0);

    if (a < 8.) {
	if (a < 1.) {
/* ----------------------------------------------------------------------- */
//                    		A < 1
/* ----------------------------------------------------------------------- */
	    if (b < 8.)
		return gamln(a) + (gamln(b) - gamln(a+b));
	    else
		return gamln(a) + algdiv(a, b);
	}
	/* else */
/* ----------------------------------------------------------------------- */
//				1 <= A < 8
/* ----------------------------------------------------------------------- */
	double w;
	if (a < 2.) {
	    if (b <= 2.) {
		return gamln(a) + gamln(b) - gsumln(a, b);
	    }
	    /* else */

	    if (b < 8.) {
		w = 0.;
		goto L40;
	    }
	    return gamln(a) + algdiv(a, b);
	}
	// else L30:    REDUCTION OF A WHEN B <= 1000

	if (b <= 1e3) {
	    int n = (int)(a - 1.);
	    w = 1.;
	    for (int i = 1; i <= n; ++i) {
		a += -1.;
		double h = a / b;
		w *= h / (h + 1.);
	    }
	    w = log(w);

	    if (b >= 8.)
		return w + gamln(a) + algdiv(a, b);

	    // else
	L40:
	    // 	1 < A <= B < 8 :  reduction of B
	    n = (int)(b - 1.);
	    double z = 1.;
	    for (int i = 1; i <= n; ++i) {
		b += -1.;
		z *= b / (a + b);
	    }
	    return w + log(z) + (gamln(a) + (gamln(b) - gsumln(a, b)));
	}
	else { // L50:	reduction of A when  B > 1000
	    int n = (int)(a - 1.);
	    w = 1.;
	    for (int i = 1; i <= n; ++i) {
		a += -1.;
		w *= a / (a / b + 1.);
	    }
	    return log(w) - n * log(b) + (gamln(a) + algdiv(a, b));
	}

    } else {
/* ----------------------------------------------------------------------- */
	// L60:			A >= 8
/* ----------------------------------------------------------------------- */

	double
	    w = bcorr(a, b),
	    h = a / b,
	    u = -(a - 0.5) * log(h / (h + 1.)),
	    v = b * alnrel(h);
	if (u > v)
	    return log(b) * -0.5 + e + w - v - u;
	else
	    return log(b) * -0.5 + e + w - u - v;
    }

} /* betaln */

static double gsumln(double a, double b)
{
/* ----------------------------------------------------------------------- */
/*          EVALUATION OF THE FUNCTION LN(GAMMA(A + B)) */
/*          FOR 1 <= A <= 2  AND  1 <= B <= 2 */
/* ----------------------------------------------------------------------- */

    double x = a + b - 2.;/* in [0, 2] */

    if (x <= 0.25)
	return gamln1(x + 1.);

    /* else */
    if (x <= 1.25)
	return gamln1(x) + alnrel(x);
    /* else x > 1.25 : */
    return gamln1(x - 1.) + log(x * (x + 1.));

} /* gsumln */

static double bcorr(double a0, double b0)
{
/* ----------------------------------------------------------------------- */

/*     EVALUATION OF  DEL(A0) + DEL(B0) - DEL(A0 + B0)  WHERE */
/*     LN(GAMMA(A)) = (A - 0.5)*LN(A) - A + 0.5*LN(2*PI) + DEL(A). */
/*     IT IS ASSUMED THAT A0 >= 8 AND B0 >= 8. */

/* ----------------------------------------------------------------------- */
    /* Initialized data */

    static double c0 = .0833333333333333;
    static double c1 = -.00277777777760991;
    static double c2 = 7.9365066682539e-4;
    static double c3 = -5.9520293135187e-4;
    static double c4 = 8.37308034031215e-4;
    static double c5 = -.00165322962780713;

    /* System generated locals */
    double ret_val, r1;

    /* Local variables */
    double a, b, c, h, t, w, x, s3, s5, x2, s7, s9, s11;
/* ------------------------ */
    a = min(a0, b0);
    b = max(a0, b0);

    h = a / b;
    c = h / (h + 1.);
    x = 1. / (h + 1.);
    x2 = x * x;

/*                SET SN = (1 - X^N)/(1 - X) */

    s3 = x + x2 + 1.;
    s5 = x + x2 * s3 + 1.;
    s7 = x + x2 * s5 + 1.;
    s9 = x + x2 * s7 + 1.;
    s11 = x + x2 * s9 + 1.;

/*                SET W = DEL(B) - DEL(A + B) */

/* Computing 2nd power */
    r1 = 1. / b;
    t = r1 * r1;
    w = ((((c5 * s11 * t + c4 * s9) * t + c3 * s7) * t + c2 * s5) * t + c1 *
	    s3) * t + c0;
    w *= c / b;

/*                   COMPUTE  DEL(A) + W */

/* Computing 2nd power */
    r1 = 1. / a;
    t = r1 * r1;
    ret_val = (((((c5 * t + c4) * t + c3) * t + c2) * t + c1) * t + c0) / a +
	    w;
    return ret_val;
} /* bcorr */

static double algdiv(double a, double b)
{
/* ----------------------------------------------------------------------- */

/*     COMPUTATION OF LN(GAMMA(B)/GAMMA(A+B)) WHEN B >= 8 */

/*                         -------- */

/*     IN THIS ALGORITHM, DEL(X) IS THE FUNCTION DEFINED BY */
/*     LN(GAMMA(X)) = (X - 0.5)*LN(X) - X + 0.5*LN(2*PI) + DEL(X). */

/* ----------------------------------------------------------------------- */

    /* Initialized data */

    static double c0 = .0833333333333333;
    static double c1 = -.00277777777760991;
    static double c2 = 7.9365066682539e-4;
    static double c3 = -5.9520293135187e-4;
    static double c4 = 8.37308034031215e-4;
    static double c5 = -.00165322962780713;

    double c, d, h, t, u, v, w, x, s3, s5, x2, s7, s9, s11;

/* ------------------------ */
    if (a > b) {
	h = b / a;
	c = 1. / (h + 1.);
	x = h / (h + 1.);
	d = a + (b - 0.5);
    }
    else {
	h = a / b;
	c = h / (h + 1.);
	x = 1. / (h + 1.);
	d = b + (a - 0.5);
    }

/* Set s<n> = (1 - x^n)/(1 - x) : */

    x2 = x * x;
    s3 = x + x2 + 1.;
    s5 = x + x2 * s3 + 1.;
    s7 = x + x2 * s5 + 1.;
    s9 = x + x2 * s7 + 1.;
    s11 = x + x2 * s9 + 1.;

/* w := Del(b) - Del(a + b) */

    t = 1./ (b * b);
    w = ((((c5 * s11 * t + c4 * s9) * t + c3 * s7) * t + c2 * s5) * t + c1 *
	    s3) * t + c0;
    w *= c / b;

/*                    COMBINE THE RESULTS */

    u = d * alnrel(a / b);
    v = a * (log(b) - 1.);
    if (u > v)
	return w - v - u;
    else
	return w - u - v;
} /* algdiv */

static double gamln(double a)
{
/* -----------------------------------------------------------------------
 *            Evaluation of  ln(gamma(a))  for positive a
 * ----------------------------------------------------------------------- */
/*     Written by Alfred H. Morris */
/*          Naval Surface Warfare Center */
/*          Dahlgren, Virginia */
/* ----------------------------------------------------------------------- */

    static double d = .418938533204673;/* d == 0.5*(LN(2*PI) - 1) */

    static double c0 = .0833333333333333;
    static double c1 = -.00277777777760991;
    static double c2 = 7.9365066682539e-4;
    static double c3 = -5.9520293135187e-4;
    static double c4 = 8.37308034031215e-4;
    static double c5 = -.00165322962780713;

    if (a <= 0.8)
	return gamln1(a) - log(a); /* ln(G(a+1)) - ln(a) == ln(G(a+1)/a) = ln(G(a)) */
    else if (a <= 2.25)
	return gamln1(a - 0.5 - 0.5);

    else if (a < 10.) {
	int i, n = (int)(a - 1.25);
	double t = a;
	double w = 1.;
	for (i = 1; i <= n; ++i) {
	    t += -1.;
	    w *= t;
	}
	return gamln1(t - 1.) + log(w);
    }
    else { /* a >= 10 */
	double t = 1. / (a * a);
	double w = (((((c5 * t + c4) * t + c3) * t + c2) * t + c1) * t + c0) / a;
	return d + w + (a - 0.5) * (log(a) - 1.);
    }
} /* gamln */
/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 */

#include "from-R.h"

double fmin2(double x, double y)
{
#ifdef IEEE_754
	if (ISNAN(x) || ISNAN(y))
		return x + y;
#endif
	return (x < y) ? x : y;
}
/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 2005-6 Morten Welinder <terra@gnome.org>
 *  Copyright (C) 2005-10 The R Foundation
 *  Copyright (C) 2006-2015 The R Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 *
 *  SYNOPSIS
 *
 *	#include <Rmath.h>
 *
 *	double pgamma (double x, double alph, double scale,
 *		       int lower_tail, int log_p)
 *
 *	double log1pmx	(double x)
 *	double lgamma1p (double a)
 *
 *	double logspace_add (double logx, double logy)
 *	double logspace_sub (double logx, double logy)
 *	double logspace_sum (double* logx, int n)
 *
 *
 *  DESCRIPTION
 *
 *	This function computes the distribution function for the
 *	gamma distribution with shape parameter alph and scale parameter
 *	scale.	This is also known as the incomplete gamma function.
 *	See Abramowitz and Stegun (6.5.1) for example.
 *
 *  NOTES
 *
 *	Complete redesign by Morten Welinder, originally for Gnumeric.
 *	Improvements (e.g. "while NEEDED_SCALE") by Martin Maechler
 *
 *  REFERENCES
 *
 */

#include "from-R.h"

/*----------- DEBUGGING -------------
 * make CFLAGS='-DDEBUG_p -g'
 * (cd `R-devel RHOME`/src/nmath; gcc -I. -I../../src/include -I../../../R/src/include  -DHAVE_CONFIG_H -fopenmp -DDEBUG_p -g -c ../../../R/src/nmath/pgamma.c -o pgamma.o)
 */

/* Scalefactor:= (2^32)^8 = 2^256 = 1.157921e+77 */
#define SQR(x) ((x)*(x))
static const double scalefactor = SQR(SQR(SQR(4294967296.0)));
#undef SQR

/* If |x| > |k| * M_cutoff,  then  log[ exp(-x) * k^x ]	 =~=  -x */
static const double M_cutoff = M_LN2 * DBL_MAX_EXP / DBL_EPSILON;/*=3.196577e18*/

/* Continued fraction for calculation of
 *    1/i + x/(i+d) + x^2/(i+2*d) + x^3/(i+3*d) + ... = sum_{k=0}^Inf x^k/(i+k*d)
 *
 * auxilary in log1pmx() and lgamma1p()
 */
static double
logcf (double x, double i, double d,
       double eps /* ~ relative tolerance */)
{
    double c1 = 2 * d;
    double c2 = i + d;
    double c4 = c2 + d;
    double a1 = c2;
    double b1 = i * (c2 - i * x);
    double b2 = d * d * x;
    double a2 = c4 * c2 - b2;

#if 0
    assert (i > 0);
    assert (d >= 0);
#endif

    b2 = c4 * b1 - i * b2;

    while (fabs(a2 * b1 - a1 * b2) > fabs(eps * b1 * b2)) {
	double c3 = c2*c2*x;
	c2 += d;
	c4 += d;
	a1 = c4 * a2 - c3 * a1;
	b1 = c4 * b2 - c3 * b1;

	c3 = c1 * c1 * x;
	c1 += d;
	c4 += d;
	a2 = c4 * a1 - c3 * a2;
	b2 = c4 * b1 - c3 * b2;

	if (fabs (b2) > scalefactor) {
	    a1 /= scalefactor;
	    b1 /= scalefactor;
	    a2 /= scalefactor;
	    b2 /= scalefactor;
	} else if (fabs (b2) < 1 / scalefactor) {
	    a1 *= scalefactor;
	    b1 *= scalefactor;
	    a2 *= scalefactor;
	    b2 *= scalefactor;
	}
    }

    return a2 / b2;
}

/* Accurate calculation of log(1+x)-x, particularly for small x.  */
double log1pmx (double x)
{
    static const double minLog1Value = -0.79149064;

    if (x > 1 || x < minLog1Value)
	return log1p(x) - x;
    else { /* -.791 <=  x <= 1  -- expand in  [x/(2+x)]^2 =: y :
	    * log(1+x) - x =  x/(2+x) * [ 2 * y * S(y) - x],  with
	    * ---------------------------------------------
	    * S(y) = 1/3 + y/5 + y^2/7 + ... = \sum_{k=0}^\infty  y^k / (2k + 3)
	   */
	double r = x / (2 + x), y = r * r;
	if (fabs(x) < 1e-2) {
	    static const double two = 2;
	    return r * ((((two / 9 * y + two / 7) * y + two / 5) * y +
			    two / 3) * y - x);
	} else {
	    static const double tol_logcf = 1e-14;
	    return r * (2 * y * logcf (y, 3, 2, tol_logcf) - x);
	}
    }
}


/* Compute  log(gamma(a+1))  accurately also for small a (0 < a < 0.5). */
double lgamma1p (double a)
{
    const double eulers_const =	 0.5772156649015328606065120900824024;

    /* coeffs[i] holds (zeta(i+2)-1)/(i+2) , i = 0:(N-1), N = 40 : */
    const int N = 40;
    static const double coeffs[40] = {
	0.3224670334241132182362075833230126e-0,/* = (zeta(2)-1)/2 */
	0.6735230105319809513324605383715000e-1,/* = (zeta(3)-1)/3 */
	0.2058080842778454787900092413529198e-1,
	0.7385551028673985266273097291406834e-2,
	0.2890510330741523285752988298486755e-2,
	0.1192753911703260977113935692828109e-2,
	0.5096695247430424223356548135815582e-3,
	0.2231547584535793797614188036013401e-3,
	0.9945751278180853371459589003190170e-4,
	0.4492623673813314170020750240635786e-4,
	0.2050721277567069155316650397830591e-4,
	0.9439488275268395903987425104415055e-5,
	0.4374866789907487804181793223952411e-5,
	0.2039215753801366236781900709670839e-5,
	0.9551412130407419832857179772951265e-6,
	0.4492469198764566043294290331193655e-6,
	0.2120718480555466586923135901077628e-6,
	0.1004322482396809960872083050053344e-6,
	0.4769810169363980565760193417246730e-7,
	0.2271109460894316491031998116062124e-7,
	0.1083865921489695409107491757968159e-7,
	0.5183475041970046655121248647057669e-8,
	0.2483674543802478317185008663991718e-8,
	0.1192140140586091207442548202774640e-8,
	0.5731367241678862013330194857961011e-9,
	0.2759522885124233145178149692816341e-9,
	0.1330476437424448948149715720858008e-9,
	0.6422964563838100022082448087644648e-10,
	0.3104424774732227276239215783404066e-10,
	0.1502138408075414217093301048780668e-10,
	0.7275974480239079662504549924814047e-11,
	0.3527742476575915083615072228655483e-11,
	0.1711991790559617908601084114443031e-11,
	0.8315385841420284819798357793954418e-12,
	0.4042200525289440065536008957032895e-12,
	0.1966475631096616490411045679010286e-12,
	0.9573630387838555763782200936508615e-13,
	0.4664076026428374224576492565974577e-13,
	0.2273736960065972320633279596737272e-13,
	0.1109139947083452201658320007192334e-13/* = (zeta(40+1)-1)/(40+1) */
    };

    const double c = 0.2273736845824652515226821577978691e-12;/* zeta(N+2)-1 */
    const double tol_logcf = 1e-14;
    double lgam;
    int i;

    if (fabs (a) >= 0.5)
	return lgammafn (a + 1);

    /* Abramowitz & Stegun 6.1.33 : for |x| < 2,
     * <==> log(gamma(1+x)) = -(log(1+x) - x) - gamma*x + x^2 * \sum_{n=0}^\infty c_n (-x)^n
     * where c_n := (Zeta(n+2) - 1)/(n+2)  = coeffs[n]
     *
     * Here, another convergence acceleration trick is used to compute
     * lgam(x) :=  sum_{n=0..Inf} c_n (-x)^n
     */
    lgam = c * logcf(-a / 2, N + 2, 1, tol_logcf);
    for (i = N - 1; i >= 0; i--)
	lgam = coeffs[i] - a * lgam;

    return (a * lgam - eulers_const) * a - log1pmx (a);
} /* lgamma1p */



/*
 * Compute the log of a sum from logs of terms, i.e.,
 *
 *     log (exp (logx) + exp (logy))
 *
 * without causing overflows and without throwing away large handfuls
 * of accuracy.
 */
double logspace_add (double logx, double logy)
{
    return fmax2 (logx, logy) + log1p (exp (-fabs (logx - logy)));
}


/*
 * Compute the log of a difference from logs of terms, i.e.,
 *
 *     log (exp (logx) - exp (logy))
 *
 * without causing overflows and without throwing away large handfuls
 * of accuracy.
 */
double logspace_sub (double logx, double logy)
{
    return logx + R_Log1_Exp(logy - logx);
}

/*
 * Compute the log of a sum from logs of terms, i.e.,
 *
 *     log (sum_i  exp (logx[i]) ) =
 *     log (e^M * sum_i  e^(logx[i] - M) ) =
 *     M + log( sum_i  e^(logx[i] - M)
 *
 * without causing overflows or throwing much accuracy.
 */
#ifdef HAVE_LONG_DOUBLE
# define EXP expl
# define LOG logl
#else
# define EXP exp
# define LOG log
#endif
double logspace_sum (const double* logx, int n)
{
    if(n == 0) return ML_NEGINF; // = log( sum(<empty>) )
    if(n == 1) return logx[0];
    if(n == 2) return logspace_add(logx[0], logx[1]);
    // else (n >= 3) :
    int i;
    // Mx := max_i log(x_i)
    double Mx = logx[0];
    for(i = 1; i < n; i++) if(Mx < logx[i]) Mx = logx[i];
    LDOUBLE s = (LDOUBLE) 0.;
    for(i = 0; i < n; i++) s += EXP(logx[i] - Mx);
    return Mx + (double) LOG(s);
}



/* dpois_wrap (x__1, lambda) := dpois(x__1 - 1, lambda);  where
 * dpois(k, L) := exp(-L) L^k / gamma(k+1)  {the usual Poisson probabilities}
 *
 * and  dpois*(.., give_log = TRUE) :=  log( dpois*(..) )
*/
static double
dpois_wrap (double x_plus_1, double lambda, int give_log)
{
#ifdef DEBUG_p
    REprintf (" dpois_wrap(x+1=%.14g, lambda=%.14g, log=%d)\n",
	      x_plus_1, lambda, give_log);
#endif
    if (!R_FINITE(lambda))
	return R_D__0;
    if (x_plus_1 > 1)
	return dpois_raw (x_plus_1 - 1, lambda, give_log);
    if (lambda > fabs(x_plus_1 - 1) * M_cutoff)
	return R_D_exp(-lambda - lgammafn(x_plus_1));
    else {
	double d = dpois_raw (x_plus_1, lambda, give_log);
#ifdef DEBUG_p
	REprintf ("  -> d=dpois_raw(..)=%.14g\n", d);
#endif
	return give_log
	    ? d + log (x_plus_1 / lambda)
	    : d * (x_plus_1 / lambda);
    }
}

/*
 * Abramowitz and Stegun 6.5.29 [right]
 */
static double
pgamma_smallx (double x, double alph, int lower_tail, int log_p)
{
    double sum = 0, c = alph, n = 0, term;

#ifdef DEBUG_p
    REprintf (" pg_smallx(x=%.12g, alph=%.12g): ", x, alph);
#endif

    /*
     * Relative to 6.5.29 all terms have been multiplied by alph
     * and the first, thus being 1, is omitted.
     */

    do {
	n++;
	c *= -x / n;
	term = c / (alph + n);
	sum += term;
    } while (fabs (term) > DBL_EPSILON * fabs (sum));

#ifdef DEBUG_p
    REprintf ("%5.0f terms --> conv.sum=%g;", n, sum);
#endif
    if (lower_tail) {
	double f1 = log_p ? log1p (sum) : 1 + sum;
	double f2;
	if (alph > 1) {
	    f2 = dpois_raw (alph, x, log_p);
	    f2 = log_p ? f2 + x : f2 * exp (x);
	} else if (log_p)
	    f2 = alph * log (x) - lgamma1p (alph);
	else
	    f2 = pow (x, alph) / exp (lgamma1p (alph));
#ifdef DEBUG_p
    REprintf (" (f1,f2)= (%g,%g)\n", f1,f2);
#endif
	return log_p ? f1 + f2 : f1 * f2;
    } else {
	double lf2 = alph * log (x) - lgamma1p (alph);
#ifdef DEBUG_p
	REprintf (" 1:%.14g  2:%.14g\n", alph * log (x), lgamma1p (alph));
	REprintf (" sum=%.14g  log(1+sum)=%.14g	 lf2=%.14g\n",
		  sum, log1p (sum), lf2);
#endif
	if (log_p)
	    return R_Log1_Exp (log1p (sum) + lf2);
	else {
	    double f1m1 = sum;
	    double f2m1 = expm1 (lf2);
	    return -(f1m1 + f2m1 + f1m1 * f2m1);
	}
    }
} /* pgamma_smallx() */

static double
pd_upper_series (double x, double y, int log_p)
{
    double term = x / y;
    double sum = term;

    do {
	y++;
	term *= x / y;
	sum += term;
    } while (term > sum * DBL_EPSILON);

    /* sum =  \sum_{n=1}^ oo  x^n     / (y*(y+1)*...*(y+n-1))
     *	   =  \sum_{n=0}^ oo  x^(n+1) / (y*(y+1)*...*(y+n))
     *	   =  x/y * (1 + \sum_{n=1}^oo	x^n / ((y+1)*...*(y+n)))
     *	   ~  x/y +  o(x/y)   {which happens when alph -> Inf}
     */
    return log_p ? log (sum) : sum;
}

/* Continued fraction for calculation of
 *    scaled upper-tail F_{gamma}
 *  ~=  (y / d) * [1 +  (1-y)/d +  O( ((1-y)/d)^2 ) ]
 */
static double
pd_lower_cf (double y, double d)
{
    double f= 0.0 /* -Wall */, of, f0;
    double i, c2, c3, c4,  a1, b1,  a2, b2;

#define	NEEDED_SCALE				\
	  (b2 > scalefactor) {			\
	    a1 /= scalefactor;			\
	    b1 /= scalefactor;			\
	    a2 /= scalefactor;			\
	    b2 /= scalefactor;			\
	}

#define max_it 200000

#ifdef DEBUG_p
    REprintf("pd_lower_cf(y=%.14g, d=%.14g)", y, d);
#endif
    if (y == 0) return 0;

    f0 = y/d;
    /* Needed, e.g. for  pgamma(10^c(100,295), shape= 1.1, log=TRUE): */
    if(fabs(y - 1) < fabs(d) * DBL_EPSILON) { /* includes y < d = Inf */
#ifdef DEBUG_p
	REprintf(" very small 'y' -> returning (y/d)\n");
#endif
	return (f0);
    }

    if(f0 > 1.) f0 = 1.;
    c2 = y;
    c4 = d; /* original (y,d), *not* potentially scaled ones!*/

    a1 = 0; b1 = 1;
    a2 = y; b2 = d;

    while NEEDED_SCALE

    i = 0; of = -1.; /* far away */
    while (i < max_it) {

	i++;	c2--;	c3 = i * c2;	c4 += 2;
	/* c2 = y - i,  c3 = i(y - i),  c4 = d + 2i,  for i odd */
	a1 = c4 * a2 + c3 * a1;
	b1 = c4 * b2 + c3 * b1;

	i++;	c2--;	c3 = i * c2;	c4 += 2;
	/* c2 = y - i,  c3 = i(y - i),  c4 = d + 2i,  for i even */
	a2 = c4 * a1 + c3 * a2;
	b2 = c4 * b1 + c3 * b2;

	if NEEDED_SCALE

	if (b2 != 0) {
	    f = a2 / b2;
	    /* convergence check: relative; "absolute" for very small f : */
	    if (fabs (f - of) <= DBL_EPSILON * fmax2(f0, fabs(f))) {
#ifdef DEBUG_p
		REprintf(" %g iter.\n", i);
#endif
		return f;
	    }
	    of = f;
	}
    }

    MATHLIB_WARNING(" ** NON-convergence in pgamma()'s pd_lower_cf() f= %g.\n",
		    f);
    return f;/* should not happen ... */
} /* pd_lower_cf() */
#undef NEEDED_SCALE


static double
pd_lower_series (double lambda, double y)
{
    double term = 1, sum = 0;

#ifdef DEBUG_p
    REprintf("pd_lower_series(lam=%.14g, y=%.14g) ...", lambda, y);
#endif
    while (y >= 1 && term > sum * DBL_EPSILON) {
	term *= y / lambda;
	sum += term;
	y--;
    }
    /* sum =  \sum_{n=0}^ oo  y*(y-1)*...*(y - n) / lambda^(n+1)
     *	   =  y/lambda * (1 + \sum_{n=1}^Inf  (y-1)*...*(y-n) / lambda^n)
     *	   ~  y/lambda + o(y/lambda)
     */
#ifdef DEBUG_p
    REprintf(" done: term=%g, sum=%g, y= %g\n", term, sum, y);
#endif

    if (y != floor (y)) {
	/*
	 * The series does not converge as the terms start getting
	 * bigger (besides flipping sign) for y < -lambda.
	 */
	double f;
#ifdef DEBUG_p
	REprintf(" y not int: add another term ");
#endif
	/* FIXME: in quite few cases, adding  term*f  has no effect (f too small)
	 *	  and is unnecessary e.g. for pgamma(4e12, 121.1) */
	f = pd_lower_cf (y, lambda + 1 - y);
#ifdef DEBUG_p
	REprintf("  (= %.14g) * term = %.14g to sum %g\n", f, term * f, sum);
#endif
	sum += term * f;
    }

    return sum;
} /* pd_lower_series() */

/*
 * Compute the following ratio with higher accuracy that would be had
 * from doing it directly.
 *
 *		 dnorm (x, 0, 1, FALSE)
 *	   ----------------------------------
 *	   pnorm (x, 0, 1, lower_tail, FALSE)
 *
 * Abramowitz & Stegun 26.2.12
 */
static double
dpnorm (double x, int lower_tail, double lp)
{
    /*
     * So as not to repeat a pnorm call, we expect
     *
     *	 lp == pnorm (x, 0, 1, lower_tail, TRUE)
     *
     * but use it only in the non-critical case where either x is small
     * or p==exp(lp) is close to 1.
     */

    if (x < 0) {
	x = -x;
	lower_tail = !lower_tail;
    }

    if (x > 10 && !lower_tail) {
	double term = 1 / x;
	double sum = term;
	double x2 = x * x;
	double i = 1;

	do {
	    term *= -i / x2;
	    sum += term;
	    i += 2;
	} while (fabs (term) > DBL_EPSILON * sum);

	return 1 / sum;
    } else {
	double d = dnorm (x, 0., 1., FALSE);
	return d / exp (lp);
    }
}

/*
 * Asymptotic expansion to calculate the probability that Poisson variate
 * has value <= x.
 * Various assertions about this are made (without proof) at
 * http://members.aol.com/iandjmsmith/PoissonApprox.htm
 */
static double
ppois_asymp (double x, double lambda, int lower_tail, int log_p)
{
    static const double coefs_a[8] = {
	-1e99, /* placeholder used for 1-indexing */
	2/3.,
	-4/135.,
	8/2835.,
	16/8505.,
	-8992/12629925.,
	-334144/492567075.,
	698752/1477701225.
    };

    static const double coefs_b[8] = {
	-1e99, /* placeholder */
	1/12.,
	1/288.,
	-139/51840.,
	-571/2488320.,
	163879/209018880.,
	5246819/75246796800.,
	-534703531/902961561600.
    };

    double elfb, elfb_term;
    double res12, res1_term, res1_ig, res2_term, res2_ig;
    double dfm, pt_, s2pt, f, np;
    int i;

    dfm = lambda - x;
    /* If lambda is large, the distribution is highly concentrated
       about lambda.  So representation error in x or lambda can lead
       to arbitrarily large values of pt_ and hence divergence of the
       coefficients of this approximation.
    */
    pt_ = - log1pmx (dfm / x);
    s2pt = sqrt (2 * x * pt_);
    if (dfm < 0) s2pt = -s2pt;

    res12 = 0;
    res1_ig = res1_term = sqrt (x);
    res2_ig = res2_term = s2pt;
    for (i = 1; i < 8; i++) {
	res12 += res1_ig * coefs_a[i];
	res12 += res2_ig * coefs_b[i];
	res1_term *= pt_ / i ;
	res2_term *= 2 * pt_ / (2 * i + 1);
	res1_ig = res1_ig / x + res1_term;
	res2_ig = res2_ig / x + res2_term;
    }

    elfb = x;
    elfb_term = 1;
    for (i = 1; i < 8; i++) {
	elfb += elfb_term * coefs_b[i];
	elfb_term /= x;
    }
    if (!lower_tail) elfb = -elfb;
#ifdef DEBUG_p
    REprintf ("res12 = %.14g   elfb=%.14g\n", elfb, res12);
#endif

    f = res12 / elfb;

    np = pnorm (s2pt, 0.0, 1.0, !lower_tail, log_p);

    if (log_p) {
	double n_d_over_p = dpnorm (s2pt, !lower_tail, np);
#ifdef DEBUG_p
	REprintf ("pp*_asymp(): f=%.14g	 np=e^%.14g  nd/np=%.14g  f*nd/np=%.14g\n",
		  f, np, n_d_over_p, f * n_d_over_p);
#endif
	return np + log1p (f * n_d_over_p);
    } else {
	double nd = dnorm (s2pt, 0., 1., log_p);

#ifdef DEBUG_p
	REprintf ("pp*_asymp(): f=%.14g	 np=%.14g  nd=%.14g  f*nd=%.14g\n",
		  f, np, nd, f * nd);
#endif
	return np + f * nd;
    }
} /* ppois_asymp() */


double pgamma_raw (double x, double alph, int lower_tail, int log_p)
{
/* Here, assume that  (x,alph) are not NA  &  alph > 0 . */

    double res;

#ifdef DEBUG_p
    REprintf("pgamma_raw(x=%.14g, alph=%.14g, low=%d, log=%d)\n",
	     x, alph, lower_tail, log_p);
#endif
    R_P_bounds_01(x, 0., ML_POSINF);

    if (x < 1) {
	res = pgamma_smallx (x, alph, lower_tail, log_p);
    } else if (x <= alph - 1 && x < 0.8 * (alph + 50)) {
	/* incl. large alph compared to x */
	double sum = pd_upper_series (x, alph, log_p);/* = x/alph + o(x/alph) */
	double d = dpois_wrap (alph, x, log_p);
#ifdef DEBUG_p
	REprintf(" alph 'large': sum=pd_upper*()= %.12g, d=dpois_w(*)= %.12g\n",
		 sum, d);
#endif
	if (!lower_tail)
	    res = log_p
		? R_Log1_Exp (d + sum)
		: 1 - d * sum;
	else
	    res = log_p ? sum + d : sum * d;
    } else if (alph - 1 < x && alph < 0.8 * (x + 50)) {
	/* incl. large x compared to alph */
	double sum;
	double d = dpois_wrap (alph, x, log_p);
#ifdef DEBUG_p
	REprintf(" x 'large': d=dpois_w(*)= %.14g ", d);
#endif
	if (alph < 1) {
	    if (x * DBL_EPSILON > 1 - alph)
		sum = R_D__1;
	    else {
		double f = pd_lower_cf (alph, x - (alph - 1)) * x / alph;
		/* = [alph/(x - alph+1) + o(alph/(x-alph+1))] * x/alph = 1 + o(1) */
		sum = log_p ? log (f) : f;
	    }
	} else {
	    sum = pd_lower_series (x, alph - 1);/* = (alph-1)/x + o((alph-1)/x) */
	    sum = log_p ? log1p (sum) : 1 + sum;
	}
#ifdef DEBUG_p
	REprintf(", sum= %.14g\n", sum);
#endif
	if (!lower_tail)
	    res = log_p ? sum + d : sum * d;
	else
	    res = log_p
		? R_Log1_Exp (d + sum)
		: 1 - d * sum;
    } else { /* x >= 1 and x fairly near alph. */
#ifdef DEBUG_p
	REprintf(" using ppois_asymp()\n");
#endif
	res = ppois_asymp (alph - 1, x, !lower_tail, log_p);
    }

    /*
     * We lose a fair amount of accuracy to underflow in the cases
     * where the final result is very close to DBL_MIN.	 In those
     * cases, simply redo via log space.
     */
    if (!log_p && res < DBL_MIN / DBL_EPSILON) {
	/* with(.Machine, double.xmin / double.eps) #|-> 1.002084e-292 */
#ifdef DEBUG_p
	REprintf(" very small res=%.14g; -> recompute via log\n", res);
#endif
	return exp (pgamma_raw (x, alph, lower_tail, 1));
    } else
	return res;
}


double pgamma(double x, double alph, double scale, int lower_tail, int log_p)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(alph) || ISNAN(scale))
	return x + alph + scale;
#endif
    if(alph < 0. || scale <= 0.)
	ML_ERR_return_NAN;
    x /= scale;
#ifdef IEEE_754
    if (ISNAN(x)) /* eg. original x = scale = +Inf */
	return x;
#endif
    if(alph == 0.) /* limit case; useful e.g. in pnchisq() */
	return (x <= 0) ? R_DT_0: R_DT_1; /* <= assert  pgamma(0,0) ==> 0 */
    return pgamma_raw (x, alph, lower_tail, log_p);
}
/* From: terra@gnome.org (Morten Welinder)
 * To: R-bugs@biostat.ku.dk
 * Cc: maechler@stat.math.ethz.ch
 * Subject: Re: [Rd] pgamma discontinuity (PR#7307)
 * Date: Tue, 11 Jan 2005 13:57:26 -0500 (EST)

 * this version of pgamma appears to be quite good and certainly a vast
 * improvement over current R code.  (I last looked at 2.0.1)  Apart from
 * type naming, this is what I have been using for Gnumeric 1.4.1.

 * This could be included into R as-is, but you might want to benefit from
 * making logcf, log1pmx, lgamma1p, and possibly logspace_add/logspace_sub
 * available to other parts of R.

 * MM: I've not (yet?) taken  logcf(), but the other four
 */
/*
 *  AUTHOR
 *    Catherine Loader, catherine@research.bell-labs.com.
 *    October 23, 2000.
 *
 *  Merge in to R:
 *	Copyright (C) 2000-2016 The R Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 *
 *
 * DESCRIPTION
 *
 *    dpois() checks argument validity and calls dpois_raw().
 *
 *    dpois_raw() computes the Poisson probability  lb^x exp(-lb) / x!.
 *      This does not check that x is an integer, since dgamma() may
 *      call this with a fractional x argument. Any necessary argument
 *      checks should be done in the calling function.
 *
 */

#include "from-R.h"


// called also from dgamma.c, pgamma.c, dnbeta.c, dnbinom.c, dnchisq.c :
double dpois_raw(double x, double lambda, int give_log)
{
    /*       x >= 0 ; integer for dpois(), but not e.g. for pgamma()!
        lambda >= 0
    */
    if (lambda == 0) return( (x == 0) ? R_D__1 : R_D__0 );
    if (!R_FINITE(lambda)) return R_D__0; // including for the case where  x = lambda = +Inf
    if (x < 0) return( R_D__0 );
    if (x <= lambda * DBL_MIN) return(R_D_exp(-lambda) );
    if (lambda < x * DBL_MIN) {
	if (!R_FINITE(x)) // lambda < x = +Inf
	    return R_D__0;
	// else
	return(R_D_exp(-lambda + x*log(lambda) -lgammafn(x+1)));
    }
    return(R_D_fexp( M_2PI*x, -stirlerr(x)-bd0(x,lambda) ));
}

double dpois(double x, double lambda, int give_log)
{
#ifdef IEEE_754
    if(ISNAN(x) || ISNAN(lambda))
        return x + lambda;
#endif

    if (lambda < 0) ML_ERR_return_NAN;
    R_D_nonint_check(x);
    if (x < 0 || !R_FINITE(x))
	return R_D__0;

    x = R_forceint(x);

    return( dpois_raw(x,lambda,give_log) );
}
/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998 Ross Ihaka
 *  Copyright (C) 2000-2014 The R Core Team
 *  Copyright (C) 2003	    The R Foundation
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 *
 *  SYNOPSIS
 *
 *	double dnorm4(double x, double mu, double sigma, int give_log)
 *	      {dnorm (..) is synonymous and preferred inside R}
 *
 *  DESCRIPTION
 *
 *	Compute the density of the normal distribution.
 */

#include "from-R.h"


double dnorm4(double x, double mu, double sigma, int give_log)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(mu) || ISNAN(sigma))
	return x + mu + sigma;
#endif
    if(!R_FINITE(sigma)) return R_D__0;
    if(!R_FINITE(x) && mu == x) return ML_NAN;/* x-mu is NaN */
    if (sigma <= 0) {
	if (sigma < 0) ML_ERR_return_NAN;
	/* sigma == 0 */
	return (x == mu) ? ML_POSINF : R_D__0;
    }
    x = (x - mu) / sigma;

    if(!R_FINITE(x)) return R_D__0;

    x = fabs (x);
    if (x >= 2 * sqrt(DBL_MAX)) return R_D__0;
    if (give_log)
	return -(M_LN_SQRT_2PI + 0.5 * x * x + log(sigma));
    //  M_1_SQRT_2PI = 1 / sqrt(2 * pi)
#ifdef MATHLIB_FAST_dnorm
    // and for R <= 3.0.x and R-devel upto 2014-01-01:
    return M_1_SQRT_2PI * exp(-0.5 * x * x) / sigma;
#else
    // more accurate, less fast :
    if (x < 5)    return M_1_SQRT_2PI * exp(-0.5 * x * x) / sigma;

    /* ELSE:

     * x*x  may lose upto about two digits accuracy for "large" x
     * Morten Welinder's proposal for PR#15620
     * https://bugs.r-project.org/bugzilla/show_bug.cgi?id=15620

     * -- 1 --  No hoop jumping when we underflow to zero anyway:

     *  -x^2/2 <         log(2)*.Machine$double.min.exp  <==>
     *     x   > sqrt(-2*log(2)*.Machine$double.min.exp) =IEEE= 37.64031
     * but "thanks" to denormalized numbers, underflow happens a bit later,
     *  effective.D.MIN.EXP <- with(.Machine, double.min.exp + double.ulp.digits)
     * for IEEE, DBL_MIN_EXP is -1022 but "effective" is -1074
     * ==> boundary = sqrt(-2*log(2)*(.Machine$double.min.exp + .Machine$double.ulp.digits))
     *              =IEEE=  38.58601
     * [on one x86_64 platform, effective boundary a bit lower: 38.56804]
     */
    if (x > sqrt(-2*M_LN2*(DBL_MIN_EXP + 1-DBL_MANT_DIG))) return 0.;

    /* Now, to get full accurary, split x into two parts,
     *  x = x1+x2, such that |x2| <= 2^-16.
     * Assuming that we are using IEEE doubles, that means that
     * x1*x1 is error free for x<1024 (but we have x < 38.6 anyway).

     * If we do not have IEEE this is still an improvement over the naive formula.
     */
    double x1 = //  R_forceint(x * 65536) / 65536 =
	ldexp( R_forceint(ldexp(x, 16)), -16);
    double x2 = x - x1;
    return M_1_SQRT_2PI / sigma *
	(exp(-0.5 * x1 * x1) * exp( (-0.5 * x2 - x1) * x2 ) );
#endif
}
/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998	    Ross Ihaka
 *  Copyright (C) 2000-2013 The R Core Team
 *  Copyright (C) 2003	    The R Foundation
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 *
 *  SYNOPSIS
 *
 *   #include <Rmath.h>
 *
 *   double pnorm5(double x, double mu, double sigma, int lower_tail,int log_p);
 *	   {pnorm (..) is synonymous and preferred inside R}
 *
 *   void   pnorm_both(double x, double *cum, double *ccum,
 *		       int i_tail, int log_p);
 *
 *  DESCRIPTION
 *
 *	The main computation evaluates near-minimax approximations derived
 *	from those in "Rational Chebyshev approximations for the error
 *	function" by W. J. Cody, Math. Comp., 1969, 631-637.  This
 *	transportable program uses rational functions that theoretically
 *	approximate the normal distribution function to at least 18
 *	significant decimal digits.  The accuracy achieved depends on the
 *	arithmetic system, the compiler, the intrinsic functions, and
 *	proper selection of the machine-dependent constants.
 *
 *  REFERENCE
 *
 *	Cody, W. D. (1993).
 *	ALGORITHM 715: SPECFUN - A Portable FORTRAN Package of
 *	Special Function Routines and Test Drivers".
 *	ACM Transactions on Mathematical Software. 19, 22-32.
 *
 *  EXTENSIONS
 *
 *  The "_both" , lower, upper, and log_p  variants were added by
 *  Martin Maechler, Jan.2000;
 *  as well as log1p() and similar improvements later on.
 *
 *  James M. Rath contributed bug report PR#699 and patches correcting SIXTEN
 *  and if() clauses {with a bug: "|| instead of &&" -> PR #2883) more in line
 *  with the original Cody code.
 */

#include "from-R.h"

double pnorm5(double x, double mu, double sigma, int lower_tail, int log_p)
{
    double p, cp;

    /* Note: The structure of these checks has been carefully thought through.
     * For example, if x == mu and sigma == 0, we get the correct answer 1.
     */
#ifdef IEEE_754
    if(ISNAN(x) || ISNAN(mu) || ISNAN(sigma))
	return x + mu + sigma;
#endif
    if(!R_FINITE(x) && mu == x) return ML_NAN;/* x-mu is NaN */
    if (sigma <= 0) {
	if(sigma < 0) ML_ERR_return_NAN;
	/* sigma = 0 : */
	return (x < mu) ? R_DT_0 : R_DT_1;
    }
    p = (x - mu) / sigma;
    if(!R_FINITE(p))
	return (x < mu) ? R_DT_0 : R_DT_1;
    x = p;

    pnorm_both(x, &p, &cp, (lower_tail ? 0 : 1), log_p);

    return(lower_tail ? p : cp);
}

#define SIXTEN	16 /* Cutoff allowing exact "*" and "/" */

void pnorm_both(double x, double *cum, double *ccum, int i_tail, int log_p)
{
/* i_tail in {0,1,2} means: "lower", "upper", or "both" :
   if(lower) return  *cum := P[X <= x]
   if(upper) return *ccum := P[X >  x] = 1 - P[X <= x]
*/
    const static double a[5] = {
	2.2352520354606839287,
	161.02823106855587881,
	1067.6894854603709582,
	18154.981253343561249,
	0.065682337918207449113
    };
    const static double b[4] = {
	47.20258190468824187,
	976.09855173777669322,
	10260.932208618978205,
	45507.789335026729956
    };
    const static double c[9] = {
	0.39894151208813466764,
	8.8831497943883759412,
	93.506656132177855979,
	597.27027639480026226,
	2494.5375852903726711,
	6848.1904505362823326,
	11602.651437647350124,
	9842.7148383839780218,
	1.0765576773720192317e-8
    };
    const static double d[8] = {
	22.266688044328115691,
	235.38790178262499861,
	1519.377599407554805,
	6485.558298266760755,
	18615.571640885098091,
	34900.952721145977266,
	38912.003286093271411,
	19685.429676859990727
    };
    const static double p[6] = {
	0.21589853405795699,
	0.1274011611602473639,
	0.022235277870649807,
	0.001421619193227893466,
	2.9112874951168792e-5,
	0.02307344176494017303
    };
    const static double q[5] = {
	1.28426009614491121,
	0.468238212480865118,
	0.0659881378689285515,
	0.00378239633202758244,
	7.29751555083966205e-5
    };

    double xden, xnum, temp, del, eps, xsq, y;
#ifdef NO_DENORMS
    double min = DBL_MIN;
#endif
    int i, lower, upper;

#ifdef IEEE_754
    if(ISNAN(x)) { *cum = *ccum = x; return; }
#endif

    /* Consider changing these : */
    eps = DBL_EPSILON * 0.5;

    /* i_tail in {0,1,2} =^= {lower, upper, both} */
    lower = i_tail != 1;
    upper = i_tail != 0;

    y = fabs(x);
    if (y <= 0.67448975) { /* qnorm(3/4) = .6744.... -- earlier had 0.66291 */
	if (y > eps) {
	    xsq = x * x;
	    xnum = a[4] * xsq;
	    xden = xsq;
	    for (i = 0; i < 3; ++i) {
		xnum = (xnum + a[i]) * xsq;
		xden = (xden + b[i]) * xsq;
	    }
	} else xnum = xden = 0.0;

	temp = x * (xnum + a[3]) / (xden + b[3]);
	if(lower)  *cum = 0.5 + temp;
	if(upper) *ccum = 0.5 - temp;
	if(log_p) {
	    if(lower)  *cum = log(*cum);
	    if(upper) *ccum = log(*ccum);
	}
    }
    else if (y <= M_SQRT_32) {

	/* Evaluate pnorm for 0.674.. = qnorm(3/4) < |x| <= sqrt(32) ~= 5.657 */

	xnum = c[8] * y;
	xden = y;
	for (i = 0; i < 7; ++i) {
	    xnum = (xnum + c[i]) * y;
	    xden = (xden + d[i]) * y;
	}
	temp = (xnum + c[7]) / (xden + d[7]);

#define do_del(X)							\
	xsq = trunc(X * SIXTEN) / SIXTEN;				\
	del = (X - xsq) * (X + xsq);					\
	if(log_p) {							\
	    *cum = (-xsq * xsq * 0.5) + (-del * 0.5) + log(temp);	\
	    if((lower && x > 0.) || (upper && x <= 0.))			\
		  *ccum = log1p(-exp(-xsq * xsq * 0.5) *		\
				exp(-del * 0.5) * temp);		\
	}								\
	else {								\
	    *cum = exp(-xsq * xsq * 0.5) * exp(-del * 0.5) * temp;	\
	    *ccum = 1.0 - *cum;						\
	}

#define swap_tail						\
	if (x > 0.) {/* swap  ccum <--> cum */			\
	    temp = *cum; if(lower) *cum = *ccum; *ccum = temp;	\
	}

	do_del(y);
	swap_tail;
    }

/* else	  |x| > sqrt(32) = 5.657 :
 * the next two case differentiations were really for lower=T, log=F
 * Particularly	 *not*	for  log_p !

 * Cody had (-37.5193 < x  &&  x < 8.2924) ; R originally had y < 50
 *
 * Note that we do want symmetry(0), lower/upper -> hence use y
 */
    else if((log_p && y < 1e170) /* avoid underflow below */
	/*  ^^^^^ MM FIXME: can speedup for log_p and much larger |x| !
	 * Then, make use of  Abramowitz & Stegun, 26.2.13, something like

	 xsq = x*x;

	 if(xsq * DBL_EPSILON < 1.)
	    del = (1. - (1. - 5./(xsq+6.)) / (xsq+4.)) / (xsq+2.);
	 else
	    del = 0.;
	 *cum = -.5*xsq - M_LN_SQRT_2PI - log(x) + log1p(-del);
	 *ccum = log1p(-exp(*cum)); /.* ~ log(1) = 0 *./

 	 swap_tail;

	 [Yes, but xsq might be infinite.]

	*/
	    || (lower && -37.5193 < x  &&  x < 8.2924)
	    || (upper && -8.2924  < x  &&  x < 37.5193)
	) {

	/* Evaluate pnorm for x in (-37.5, -5.657) union (5.657, 37.5) */
	xsq = 1.0 / (x * x); /* (1./x)*(1./x) might be better */
	xnum = p[5] * xsq;
	xden = xsq;
	for (i = 0; i < 4; ++i) {
	    xnum = (xnum + p[i]) * xsq;
	    xden = (xden + q[i]) * xsq;
	}
	temp = xsq * (xnum + p[4]) / (xden + q[4]);
	temp = (M_1_SQRT_2PI - temp) / y;

	do_del(x);
	swap_tail;
    } else { /* large x such that probs are 0 or 1 */
	if(x > 0) {	*cum = R_D__1; *ccum = R_D__0;	}
	else {	        *cum = R_D__0; *ccum = R_D__1;	}
    }


#ifdef NO_DENORMS
    /* do not return "denormalized" -- we do in R */
    if(log_p) {
	if(*cum > -min)	 *cum = -0.;
	if(*ccum > -min)*ccum = -0.;
    }
    else {
	if(*cum < min)	 *cum = 0.;
	if(*ccum < min)	*ccum = 0.;
    }
#endif
    return;
}
/*
 *  Mathlib - A Mathematical Function Library
 *  Copyright (C) 1998  Ross Ihaka
 *  Copyright (C) 2000-7 The R Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 */

#include "from-R.h"
#include <limits.h>

attribute_hidden int Rf_i1mach(int i)
{
    switch(i) {

    case  1: return 5;
    case  2: return 6;
    case  3: return 0;
    case  4: return 0;

    case  5: return CHAR_BIT * sizeof(int);
    case  6: return sizeof(int)/sizeof(char);

    case  7: return 2;
    case  8: return CHAR_BIT * sizeof(int) - 1;
    case  9: return INT_MAX;

    case 10: return FLT_RADIX;

    case 11: return FLT_MANT_DIG;
    case 12: return FLT_MIN_EXP;
    case 13: return FLT_MAX_EXP;

    case 14: return DBL_MANT_DIG;
    case 15: return DBL_MIN_EXP;
    case 16: return DBL_MAX_EXP;

    default: return 0;
    }
}

int F77_NAME(i1mach)(int *i)
{
    return Rf_i1mach(*i);
}
/*
 *  Mathlib - A Mathematical Function Library
 *  Copyright (C) 1998  Ross Ihaka
 *  Copyright (C) 2000-2014 The R Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 */

/* NaNs propagated correctly */


/*-- FIXME:  Eliminate calls to these
 *   =====   o   from C code when
 *	     o   it is only used to initialize "static" variables (threading)
 *  and use the DBL_... constants instead
 */

#include "from-R.h"

attribute_hidden double Rf_d1mach(int i)
{
    switch(i) {
    case 1: return DBL_MIN;
    case 2: return DBL_MAX;

    case 3: /* = FLT_RADIX  ^ - DBL_MANT_DIG
	      for IEEE:  = 2^-53 = 1.110223e-16 = .5*DBL_EPSILON */
	return 0.5*DBL_EPSILON;

    case 4: /* = FLT_RADIX  ^ (1- DBL_MANT_DIG) =
	      for IEEE:  = 2^-52 = DBL_EPSILON */
	return DBL_EPSILON;

    case 5: return M_LOG10_2;

    default: return 0.0;
    }
}

#ifdef __cplusplus
extern "C" 
#endif

double F77_NAME(d1mach)(int *i)
{
    return Rf_d1mach(*i);
}
#include "from-R.h"
int R_finite(double x)
{
	return R_FINITE(x);
}

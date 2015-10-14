#if defined(INLA_EXPERIMENTAL)

/* R-interface.c
 * 
 * Copyright (C) 2014 Havard Rue
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
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
 *       H{\aa}vard Rue
 *       Department of Mathematical Sciences
 *       The Norwegian University of Science and Technology
 *       N-7491 Trondheim, Norway
 *       Voice: +47-7359-3533    URL  : http://www.math.ntnu.no/~hrue  
 *       Fax  : +47-7359-3524    Email: havard.rue@math.ntnu.no
 *
 */
#ifndef HGVERSION
#define HGVERSION
#endif
static const char RCSId[] = HGVERSION;

#include <stddef.h>
#include <stdio.h>
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif
#include <stdlib.h>
#include <omp.h>

#include <R.h>
#include <Rembedded.h>
#include <Rinternals.h>
#include <Rdefines.h>
#define CSTACK_DEFNS 1
#include <Rinterface.h>

//#include "GMRFLib/GMRFLib.h"
//#include "GMRFLib/GMRFLibP.h"

//#include "inla.h"
#include "R-interface.h"

#define INLA_OK (1)

static int R_init = !INLA_OK;
static int R_debug = !INLA_OK;

void inla_R_exit(void)
{
	if (R_debug)
		fprintf(stderr, "R-interface: exit\n");
#pragma omp critical
	{
		Rf_endEmbeddedR(0);
		R_init = !INLA_OK;
	}
}
int inla_R_init(void)
{
	if (R_init == !INLA_OK) {
#pragma omp critical
		{
			if (R_init == !INLA_OK) {
				char *Rargv[] = { "REmbeddedPostgres", "--gui=none", "--silent", "--vanilla" };
				int Rargc = sizeof(Rargv) / sizeof(Rargv[0]);
				Rf_initEmbeddedR(Rargc, Rargv);

				// Disable C stack limit check
				R_CStackLimit = (uintptr_t) - 1;
				R_init = INLA_OK;
				if (R_debug)
					fprintf(stderr, "R-interface: init\n");
			}
		}
	}

	return INLA_OK;
}
int inla_R_source(const char *filename)
{
	if (!filename)
		return INLA_OK;
	inla_R_init();

#pragma omp critical
	{
		SEXP e, result;
		int error;

		if (R_debug)
			fprintf(stderr, "R-interface: source file [%s]\n", filename);

		PROTECT(e = lang2(install("source"), mkString(filename)));
		PROTECT(result = R_tryEval(e, R_GlobalEnv, &error));
		if (error) {
			fprintf(stderr, "\n *** ERROR ***: source R-file [%s] failed.\n", filename);
			exit(1);
		}
		UNPROTECT(2);

	}
	return INLA_OK;
}
int inla_R_funcall1(int *n_out, double **x_out, const char *function, int n, double *x)
{
	return inla_R_funcall2(n_out, x_out, function, NULL, n, x);
}
int inla_R_funcall2(int *n_out, double **x_out, const char *function, const char *tag, int n, double *x)
{
	/*
	 * Call function(tag,x), where x is a double vector of length n. output is 'x_out' with length 'n_out'
	 */

	inla_R_init();
#pragma omp critical
	{
		if (R_debug)
			fprintf(stderr, "R-interface[%1d]: funcall2: function [%s] tag [%s] n [%1d]\n", omp_get_thread_num(), function, tag, n);

		int error, i;
		SEXP yy, xx, result, e;

		PROTECT(yy = mkString((tag ? tag : "<<<NoTag>>>")));
		PROTECT(xx = allocVector(REALSXP, n));
		for (i = 0; i < n; i++) {
			REAL(xx)[i] = x[i];
		}
		if (tag) {
			PROTECT(e = lang3(install(function), yy, xx));
		} else {
			PROTECT(e = lang2(install(function), xx));
		}
		PROTECT(result = R_tryEval(e, R_GlobalEnv, &error));
		if (error) {
			fprintf(stderr, "\n *** ERROR *** Calling R-function [%s] with tag [%s] and [%1d] arguments\n", function, tag, n);
			exit(1);
		}
		*n_out = (int) XLENGTH(result);
		*x_out = (double *) calloc((size_t) * n_out, sizeof(double));	/* otherwise I'' use the R-version... */
		for (i = 0; i < *n_out; i++) {
			(*x_out)[i] = REAL(result)[i];
		}
		UNPROTECT(4);
	}
	return INLA_OK;
}

#else

#include <stdlib.h>
#include <stdio.h>
#define ERROR_MESSAGE \
	if (1) {							\
		fprintf(stderr, "\n\n\n\t*** ERROR ***: The R-interface was not enabled at compile time: please recompile.\n\n\n"); \
		abort();						\
		exit(1);						\
	}


void inla_R_exit(void)
{
	ERROR_MESSAGE;
}

int inla_R_init(void)
{
	ERROR_MESSAGE;
}

int inla_R_source(const char *filename)
{
	ERROR_MESSAGE;
}

int inla_R_funcall1(int *n_out, double **x_out, const char *function, int n, double *x)
{
	ERROR_MESSAGE;
}

int inla_R_funcall2(int *n_out, double **x_out, const char *function, const char *tag, int n, double *x)
{
	ERROR_MESSAGE;
}
#undef ERROR_MESSAGE

#endif

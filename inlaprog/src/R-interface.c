#if defined(INLA_EXPERIMENTAL)

/* R-interface.c
 * 
 * Copyright (C) 2014-2016 Havard Rue
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

#include <assert.h>
#include <time.h>
#include <stdint.h>
#include <stddef.h>
#include <stdio.h>
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif
#include <stdlib.h>
#include <omp.h>

#define CSTACK_DEFNS 1
#include <R.h>
#include <Rembedded.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <Rinterface.h>

// two copies...
#define R_GENERIC_WRAPPER "inla.rgeneric.wrapper"
#define INLA_OK (0)
int my_file_exists(const char *filename);
int my_setenv(char *str, int prefix);
int GMRFLib_sprintf(char **ptr, const char *fmt, ...);

#include "R-interface.h"
static int R_init = 1;
static int R_debug = 0;

void inla_R_exit(void)
{
	if (R_debug)
		fprintf(stderr, "R-interface: exit\n");
	Rf_endEmbeddedR(0);
	R_init = 1;
}

int inla_R_init(void)
{
	if (R_init) {
		/*
		 * Check if R_HOME is set. If not, try to guess it, otherwise fail.
		 */
		char *rhome = getenv((const char *) "R_HOME");
		if (!rhome) {
			if (my_file_exists("/Library/Frameworks/R.framework/Resources") == INLA_OK) {
				GMRFLib_sprintf(&rhome, "R_HOME=/Library/Frameworks/R.framework/Resources");
			} else if (my_file_exists("/usr/lib64/R") == INLA_OK) {
				GMRFLib_sprintf(&rhome, "R_HOME=/usr/lib64/R");
			} else if (my_file_exists("/usr/lib/R") == INLA_OK) {
				GMRFLib_sprintf(&rhome, "R_HOME=/usr/lib/R");
			} else if (my_file_exists("/usr/lib32/R") == INLA_OK) {
				GMRFLib_sprintf(&rhome, "R_HOME=/usr/lib32/R");
			} else {
				fprintf(stderr, "\n\n");
				fprintf(stderr, "*** R-interface  ERROR: Environment variable R_HOME is not set.\n");
				fprintf(stderr, "*** R_interface  ERROR: Evaluate this in R:  Sys.getenv(\"R_HOME\")\n");
				fprintf(stderr, "\n\n");
				fflush(stderr);
				exit(1);
			}
			fprintf(stderr, "\n\n");
			fprintf(stderr, "*** R-interface WARNING: Environment variable R_HOME is not set.\n");
			fprintf(stderr, "*** R-interface WARNING: Set it to a _GUESSED_ value [%s]\n\n", rhome);
			fflush(stderr);
			my_setenv(rhome, 0);
			Free(rhome);
		}
		// char *Rargv[] = { "REmbeddedPostgres", "--gui=none", "--silent", "--no-save" };
		char *Rargv[] = { "REmbeddedPostgres", "--gui=none", "--no-save" };
		int Rargc = sizeof(Rargv) / sizeof(Rargv[0]);
		Rf_initEmbeddedR(Rargc, Rargv);

		// Disable C stack limit check
		R_CStackLimit = (uintptr_t) (-1);
		R_init = 0;
		if (R_debug)
			fprintf(stderr, "R-interface: init\n");
	}

	return (INLA_OK);
}

int inla_R_library(const char *library)
{
	if (!library)
		return (INLA_OK);
	inla_R_init();

	SEXP e, result, yy;
	int error;

	if (R_debug)
		fprintf(stderr, "R-interface: load library [%s]\n", library);
	yy = PROTECT(mkString(library));
	e = PROTECT(lang2(install("library"), yy));
	result = PROTECT(R_tryEval(e, R_GlobalEnv, &error));
	if (error) {
		fprintf(stderr, "\n *** ERROR ***: load library [%s] failed.\n", library);
		exit(1);
	}
	UNPROTECT(3);

	return (INLA_OK);
}

int inla_R_source(const char *filename)
{
	if (!filename)
		return (INLA_OK);
	inla_R_init();

	SEXP e, result, yy, false, true;
	int error;

	if (R_debug)
		fprintf(stderr, "R-interface: source file [%s]\n", filename);
	false = PROTECT(ScalarLogical(FALSE));
	true = PROTECT(ScalarLogical(TRUE));
	yy = PROTECT(mkString(filename));
	e = PROTECT(lang4(install("source"), yy, false, true));
	result = PROTECT(R_tryEval(e, R_GlobalEnv, &error));
	if (error) {
		fprintf(stderr, "\n *** ERROR ***: source R-file [%s] failed.\n", filename);
		exit(1);
	}
	UNPROTECT(5);


	return (INLA_OK);
}
int inla_R_load(const char *filename)
{
	if (!filename)
		return (INLA_OK);
	inla_R_init();

	SEXP e, result, yy;
	int error;

	if (R_debug)
		fprintf(stderr, "R-interface: loading file [%s]\n", filename);

	yy = PROTECT(mkString(filename));
	e = PROTECT(lang2(install("load"), yy));
	result = PROTECT(R_tryEval(e, R_GlobalEnv, &error));
	if (error) {
		fprintf(stderr, "\n *** ERROR ***: loading RData-file [%s] failed.\n", filename);
		exit(1);
	}
	UNPROTECT(3);

	return (INLA_OK);
}

int inla_R_inlaload(const char *filename)
{
	if (!filename)
		return (INLA_OK);
	inla_R_init();

	SEXP e, result, yy;
	int error;

	if (R_debug)
		fprintf(stderr, "R-interface: inla.load file [%s]\n", filename);

	yy = PROTECT(mkString(filename));
	e = PROTECT(lang2(install("inla.load"), yy));
	result = PROTECT(R_tryEval(e, R_GlobalEnv, &error));
	if (error) {
		fprintf(stderr, "\n *** ERROR ***: inla.load file [%s] failed.\n", filename);
		exit(1);
	}
	UNPROTECT(3);

	return (INLA_OK);
}

int inla_R_funcall2(int *n_out, double **x_out, const char *function, const char *tag, int n, double *x)
{
	/*
	 * Call function(tag,x), where x is a double vector of length n. output is 'x_out' with length 'n_out'
	 */

	inla_R_init();

	if (R_debug)
		fprintf(stderr, "R-interface[%1d]: funcall2: function [%s] tag [%s] n [%1d]\n", omp_get_thread_num(), function, tag, n);

	int error, i;
	SEXP xx, result, e;

	xx = PROTECT(allocVector(REALSXP, n));
	for (i = 0; i < n; i++) {
		REAL(xx)[i] = x[i];
	}
	if (tag) {
		SEXP yy = PROTECT(mkString(tag));
		e = PROTECT(lang3(install(function), yy, xx));
	} else {
		e = PROTECT(lang2(install(function), xx));
	}
	result = PROTECT(R_tryEval(e, R_GlobalEnv, &error));
	if (error) {
		fprintf(stderr, "\n *** ERROR *** Calling R-function [%s] with tag [%s] and [%1d] arguments\n", function, tag, n);
		exit(1);
	}
	*n_out = (int) XLENGTH(result);
	assert(*n_out >= 0);
	if (*n_out > 0) {
		*x_out = (double *) calloc((size_t) (*n_out), sizeof(double));	/* otherwise I'' use the R-version... */
		for (i = 0; i < *n_out; i++) {
			(*x_out)[i] = REAL(result)[i];
		}
	} else {
		*x_out = NULL;
	}

	UNPROTECT((tag ? 4 : 3));
	return (INLA_OK);
}

int inla_R_funcall1(int *n_out, double **x_out, const char *function, int n, double *x)
{
	return inla_R_funcall2(n_out, x_out, function, NULL, n, x);
}

int inla_R_assign(const char *variable, int n, double *x)
{
	/*
	 * variable = x
	 */

	inla_R_init();
	if (R_debug)
		fprintf(stderr, "R-interface[%1d]: assign: [%s] n [%1d]\n", omp_get_thread_num(), variable, n);

	int error, i;
	SEXP xx, result, e, yy;

	yy = PROTECT(mkString(variable));
	xx = PROTECT(allocVector(REALSXP, n));
	for (i = 0; i < n; i++) {
		REAL(xx)[i] = x[i];
	}
	e = PROTECT(lang3(install("assign"), yy, xx));
	result = PROTECT(R_tryEval(e, R_GlobalEnv, &error));
	if (error) {
		fprintf(stderr, "\n *** ERROR *** assign [%s] with n [%1d] failed\n", variable, n);
		exit(1);
	}
	UNPROTECT(4);

	return (INLA_OK);
}

int inla_R_get(int *n_out, double **x_out, const char *variable)
{
	/*
	 * return variable
	 */

	inla_R_init();
	if (R_debug)
		fprintf(stderr, "R-interface[%1d]: get: [%s]\n", omp_get_thread_num(), variable);

	int error, i;
	SEXP result, e, yy;

	yy = PROTECT(mkString(variable));
	e = PROTECT(lang2(install("get"), yy));
	result = PROTECT(R_tryEval(e, R_GlobalEnv, &error));
	if (error) {
		fprintf(stderr, "\n *** ERROR *** get [%s]\n", variable);
		exit(1);
	}
	*n_out = (int) XLENGTH(result);
	assert(*n_out >= 0);
	if (*n_out > 0) {
		*x_out = (double *) calloc((size_t) (*n_out), sizeof(double));	/* otherwise I'' use the R-version... */
		for (i = 0; i < *n_out; i++) {
			(*x_out)[i] = REAL(result)[i];
		}
	} else {
		*x_out = NULL;
	}

	UNPROTECT(3);

	return (INLA_OK);
}

int inla_R_rgeneric(int *n_out, double **x_out, const char *cmd, const char *model, int n, double *theta)
{
	/*
	 * do the rgeneric call with CMD and THETA and given MODEL name for the model definition
	 */

	inla_R_init();
	if (R_debug)
		fprintf(stderr, "R-interface[%1d]: rgeneric: [%s] model [%s]\n", omp_get_thread_num(), cmd, model);

	int error, i;
	SEXP xx_theta, result, e, yy, yyy;

	PROTECT(xx_theta = allocVector(REALSXP, n));
	for (i = 0; i < n; i++) {
		REAL(xx_theta)[i] = theta[i];
	}

	yy = PROTECT(mkString(cmd));
	yyy = PROTECT(mkString(model));
	e = PROTECT(lang4(install(R_GENERIC_WRAPPER), yy, yyy, xx_theta));
	result = PROTECT(R_tryEval(e, R_GlobalEnv, &error));
	if (error) {
		fprintf(stderr, "\n *** ERROR *** rgeneric [%s] with model [%s] failed\n", cmd, model);
		exit(1);
	}
	*n_out = (int) XLENGTH(result);
	assert(*n_out >= 0);
	if (*n_out > 0) {
		*x_out = (double *) calloc((size_t) (*n_out), sizeof(double));	/* otherwise I'' use the R-version... */
		for (i = 0; i < *n_out; i++) {
			(*x_out)[i] = REAL(result)[i];
		}
	} else {
		*x_out = NULL;
	}

	UNPROTECT(5);

	return (INLA_OK);
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

int inla_R_library(const char *filename)
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

int inla_R_assign(const char *variable, int n, double *x)
{
	ERROR_MESSAGE;
}

int inla_R_get(int *n_out, double **x_out, const char *variable)
{
	ERROR_MESSAGE;
}

int inla_R_load(const char *filename)
{
	ERROR_MESSAGE;
}

int inla_R_rgeneric(int *n_out, double **x_out, const char *cmd, const char *model, int n, double *theta)
{
	ERROR_MESSAGE;
}
int inla_R_inlaload(const char *filename)
{
	ERROR_MESSAGE;
}


#undef ERROR_MESSAGE

#endif

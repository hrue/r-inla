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
static const char RCSId[] = HGVERSION;

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
		//char *Rargv[] = { "REmbeddedPostgres", "--gui=none", "--silent", "--vanilla" };
		char *Rargv[] = { "REmbeddedPostgres", "--gui=none", "--silent", "--no-save" }; 
		int Rargc = sizeof(Rargv) / sizeof(Rargv[0]);
		Rf_initEmbeddedR(Rargc, Rargv);

		// Disable C stack limit check
		R_CStackLimit = (uintptr_t) - 1;
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

	SEXP e, result;
	int error;

	if (R_debug)
		fprintf(stderr, "R-interface: load library [%s]\n", library);
	PROTECT(e = lang2(install("library"), mkString(library)));
	PROTECT(result = R_tryEval(e, R_GlobalEnv, &error));
	if (error) {
		fprintf(stderr, "\n *** ERROR ***: load library [%s] failed.\n", library);
		exit(1);
	}
	UNPROTECT(2);

	return (INLA_OK);
}
int inla_R_source(const char *filename)
{
	if (!filename)
		return (INLA_OK);
	inla_R_init();

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


	return (INLA_OK);
}
int inla_R_load(const char *filename)
{
	if (!filename)
		return (INLA_OK);
	inla_R_init();

	SEXP e, result;
	int error;

	if (R_debug)
		fprintf(stderr, "R-interface: loading file [%s]\n", filename);

	PROTECT(e = lang2(install("load"), mkString(filename)));
	PROTECT(result = R_tryEval(e, R_GlobalEnv, &error));
	if (error) {
		fprintf(stderr, "\n *** ERROR ***: loading RData-file [%s] failed.\n", filename);
		exit(1);
	}
	UNPROTECT(2);

	return (INLA_OK);
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
	*x_out = (double *) calloc((size_t)(*n_out), sizeof(double));	/* otherwise I'' use the R-version... */
	for (i = 0; i < *n_out; i++) {
		(*x_out)[i] = REAL(result)[i];
	}
	UNPROTECT(4);

	return (INLA_OK);
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
	SEXP xx, result, e;

	PROTECT(xx = allocVector(REALSXP, n));
	for (i = 0; i < n; i++) {
		REAL(xx)[i] = x[i];
	}
	PROTECT(e = lang3(install("assign"), mkString(variable), xx));
	PROTECT(result = R_tryEval(e, R_GlobalEnv, &error));
	if (error) {
		fprintf(stderr, "\n *** ERROR *** assign [%s] with n [%1d] failed\n", variable,  n);
		exit(1);
	}
	UNPROTECT(3);

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
	SEXP result, e;

	PROTECT(e = lang2(install("get"), mkString(variable)));
	PROTECT(result = R_tryEval(e, R_GlobalEnv, &error));
	if (error) {
		fprintf(stderr, "\n *** ERROR *** get [%s]\n", variable);
		exit(1);
	}
	*n_out = (int) XLENGTH(result);
	*x_out = (double *) calloc((size_t)(*n_out), sizeof(double));	/* otherwise I'' use the R-version... */
	for (i = 0; i < *n_out; i++) {
		(*x_out)[i] = REAL(result)[i];
	}
	UNPROTECT(2);

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
	SEXP xx_theta, result, e;

	PROTECT(xx_theta = allocVector(REALSXP, n));
	for (i = 0; i < n; i++) {
		REAL(xx_theta)[i] = theta[i];
	}
	PROTECT(e = lang4(install(R_GENERIC_WRAPPER), mkString(cmd), mkString(model), xx_theta));
	PROTECT(result = R_tryEval(e, R_GlobalEnv, &error));
	if (error) {
		fprintf(stderr, "\n *** ERROR *** rgeneric [%s] with model [%s] failed\n", cmd,  model);
		exit(1);
	}
	*n_out = (int) XLENGTH(result);
	*x_out = (double *) calloc((size_t)(*n_out), sizeof(double));	/* otherwise I'' use the R-version... */
	for (i = 0; i < *n_out; i++) {
		(*x_out)[i] = REAL(result)[i];
	}
	UNPROTECT(3);

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

#undef ERROR_MESSAGE

#endif

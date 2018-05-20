
/* eval.c
 * 
 * Copyright (C) 2011  Havard Rue
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
 *        Haavard Rue
 *        CEMSE Division
 *        King Abdullah University of Science and Technology
 *        Thuwal 23955-6900, Saudi Arabia
 *        Email: haavard.rue@kaust.edu.sa
 *        Office: +966 (0)12 808 0640
 *
 */
#ifndef HGVERSION
#define HGVERSION
#endif
static const char RCSId[] = "file: " __FILE__ "  " HGVERSION;

#if !defined(__FreeBSD__)
#include <malloc.h>
#endif
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <string.h>
#include <stdio.h>

#include <muParser/muParserDLL.h>

#include "GMRFLib/GMRFLib.h"
#include "inla.h"
#include "interpol.h"
#include "eval.h"

/* 
   This is the interface to the muparser-library (http://muparser.sourceforge.net).  This code is buildt upon example2.c in the muparser/samples directory.
 */

typedef struct {
	muFloat_t **value;
	char **name;
	int n;
	int n_alloc;
	double default_value;
} eval_keep_vars_tp;

static unsigned char debug = 0;

/* 
   local functions...
 */
double inla_eval_Gamma(double arg);
double inla_eval_LogGamma(double arg);
double inla_eval_Return(double v);
double inla_eval_Not(double v);
void inla_eval_OnError(muParserHandle_t hParser);
muFloat_t *inla_eval_AddVariable(const muChar_t * a_szName, void *pUserData);
double inla_eval_lgamma(double arg);
double inla_eval_digamma(double arg);
double inla_eval_trigamma(double arg);

double inla_eval_Gamma(double arg)
{
	return exp(gsl_sf_lngamma(arg));
}
double inla_eval_LogGamma(double arg)
{
	return gsl_sf_lngamma(arg);
}
double inla_eval_lgamma(double arg)
{
	return gsl_sf_lngamma(arg);
}
double inla_eval_digamma(double arg)
{
	return gsl_sf_psi(arg);
}
double inla_eval_trigamma(double arg)
{
	return gsl_sf_psi_1(arg);
}
double inla_eval_Return(double v)
{
	return v;
}

double inla_eval_Not(double v)
{
	return v == 0.0;
}

void inla_eval_OnError(muParserHandle_t hParser)
{
	fprintf(stderr, "\n\n\nEval: Error while parsing expression:\n");
	fprintf(stderr, "----------\n");
	fprintf(stderr, "Message   :  \"%s\"\n", mupGetErrorMsg(hParser));
	fprintf(stderr, "Token     :    \"%s\"\n", mupGetErrorToken(hParser));
	fprintf(stderr, "Position  : %1d\n", mupGetErrorPos(hParser));
	fprintf(stderr, "Error Code:     %1d\n", mupGetErrorCode(hParser));
	exit(1);
}

muFloat_t *inla_eval_AddVariable(const muChar_t * a_szName, void *pUserData)
{
	eval_keep_vars_tp **aa = (eval_keep_vars_tp **) pUserData;
	if (*aa == NULL) {
		*aa = Calloc(1, eval_keep_vars_tp);
	}
	eval_keep_vars_tp *a = *aa;

	if (a->n >= a->n_alloc) {
		a->n_alloc += 16;
		if (a->n_alloc == 0) {
			assert(a->name == NULL);
			assert(a->value == NULL);
		}
		a->name = Realloc(a->name, a->n_alloc, char *);
		a->value = Realloc(a->value, a->n_alloc, muFloat_t *);
	}
	a->value[a->n] = Calloc(1, muFloat_t);
	a->value[a->n][0] = a->default_value;
	a->name[a->n] = GMRFLib_strdup(a_szName);

	if (debug) {
		printf("Eval: Add variable [%s] = %g\n", a->name[a->n], a->value[a->n][0]);
	}
	a->n++;

	return &(a->value[a->n - 1][0]);
}

double inla_eval(char *expression, double *x, double *theta, int ntheta)
{
	if (debug) {
		printf("call inla_eval with %s\n", expression);
	}

	if (strncasecmp(expression, "EXPRESSION:", strlen("EXPRESSION:")) == 0) {
		return (inla_eval_expression(expression + strlen("EXPRESSION:"), x, theta, ntheta));
	} else if (strncasecmp(expression, "TABLE:", strlen("TABLE:")) == 0) {
		return (inla_eval_table(expression + strlen("TABLE:"), x, theta, ntheta));
	} else {
		assert(0 == 1);
	}
}
double inla_eval_expression(char *expression, double *x, double *theta, int ntheta)
{
	double value;
	int i;

	/*
	 * I need this until the muparser-library is thread-safe....
	 */
#pragma omp critical
	{
		muParserHandle_t hParser;

		if (debug) {
			printf("Eval: expression: %s\n", expression);
			printf("Eval: value: %g\n", *x);
		}

		hParser = mupCreate(muBASETYPE_FLOAT);
		mupSetErrorHandler(hParser, inla_eval_OnError);

		mupSetArgSep(hParser, ';');
		mupSetDecSep(hParser, '.');
		mupSetThousandsSep(hParser, 0);
		mupDefineConst(hParser, "pi", M_PI);
		mupDefineInfixOprt(hParser, "!", inla_eval_Not, 0);
		mupDefineFun1(hParser, "return", inla_eval_Return, 1);
		mupDefineFun1(hParser, "gamma", inla_eval_Gamma, 1);
		mupDefineFun1(hParser, "lgamma", inla_eval_LogGamma, 1);
		mupDefineFun1(hParser, "digamma", inla_eval_digamma, 1);
		mupDefineFun1(hParser, "trigamma", inla_eval_trigamma, 1);
		mupDefineFun1(hParser, "log", log, 1);
		mupDefineFun1(hParser, "log10", log10, 1);
		mupDefineFun1(hParser, "ln", log, 1);
		mupDefineFun2(hParser, "pow", pow, 1);

		// add constants like THETA0, THETA1, THETA2, ...
		for (i = 0; i < ntheta; i++) {
			char *var = NULL;
			GMRFLib_sprintf(&var, "THETA%1d", i);
			mupDefineConst(hParser, var, theta[i]);
			Free(var);
		}
		eval_keep_vars_tp *keep_vars = NULL;

		keep_vars = Calloc(1, eval_keep_vars_tp);
		keep_vars->default_value = *x;
		mupSetVarFactory(hParser, inla_eval_AddVariable, (void *) &keep_vars);
		mupSetExpr(hParser, (muChar_t *) expression);
		value = (double) mupEval(hParser);

		mupRelease(hParser);

		int i;
		for (i = 0; i < keep_vars->n; i++) {
			if (debug) {
				printf("Free %s = %g\n", keep_vars->name[i], (double) keep_vars->value[i][0]);
			}
			Free(keep_vars->value[i]);
			Free(keep_vars->name[i]);
		}
		Free(keep_vars->value);
		Free(keep_vars->name);
		Free(keep_vars);
	}

	if (ISINF(value) || ISNAN(value)) {
		char *msg;

		GMRFLib_sprintf(&msg, "Expression[%s] evaluate to INF or NAN [%g] for x=%g\n", expression, value, *x);
		inla_error_general(msg);
	}

	return value;
}
double inla_eval_table(char *expression, double *xval, double *theta, int ntheta)
{
	double value;
	GMRFLib_spline_tp *s;
	GMRFLib_matrix_tp *M = NULL;

	while (*expression == ' ' || *expression == '\t') {
		expression++;
	}
	if (debug) {
		fprintf(stderr, "OPEN FILE[%s]\n", expression);
	}
	M = GMRFLib_read_fmesher_file((const char *) expression, 0, -1);
	assert(M->nrow >= 4);
	assert(M->ncol == 2);

	s = inla_spline_create(M->A, M->A + M->nrow, M->nrow);
	value = inla_spline_eval(*xval, s);

	if (0) {
		// a check of the interpolation
#pragma omp critical
		{
			double xx;
			for (xx = -20; xx < 20; xx += .1)
				printf("TABLE %g %g\n", xx, inla_spline_eval(xx, s));
			exit(1);
		}
	}

	if (ISNAN(value)) {
		char *msg;
		GMRFLib_sprintf(&msg, "table-prior returns NAN. Argument is %g but prior is defined on [%g,%g] only.", *xval, s->xmin, s->xmax);
		inla_error_general(msg);
		exit(1);
	}

	GMRFLib_matrix_free(M);
	inla_spline_free(s);

	return value;
}

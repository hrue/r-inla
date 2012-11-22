
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

double inla_eval(char *expression, double *x)
{
	if (debug) {
		printf("call inla_eval with %s\n", expression);
	}

	if (strncasecmp(expression, "EXPRESSION:", strlen("EXPRESSION:")) == 0) {
		return (inla_eval_expression(expression + strlen("EXPRESSION:"), x));
	} else if (strncasecmp(expression, "TABLE:", strlen("TABLE:")) == 0) {
		return (inla_eval_table(expression + strlen("TABLE:"), x));
	} else {
		assert(0 == 1);
	}
}
double inla_eval_expression(char *expression, double *x)
{
	double value;

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

		hParser = mupCreate();
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

		eval_keep_vars_tp *keep_vars;

		keep_vars = Calloc(1, eval_keep_vars_tp);
		keep_vars->default_value = *x;
		mupSetVarFactory(hParser, inla_eval_AddVariable, (void *) &keep_vars);
		mupSetExpr(hParser, (muChar_t *) expression);
		value = (double) mupEval(hParser);

		mupRelease(hParser);
		if (keep_vars) {
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
	}

	if (ISINF(value) || ISNAN(value)) {
		char *msg;

		GMRFLib_sprintf(&msg, "Expression[%s] evaluate to INF or NAN [%g] for x=%g\n", expression, value, *x);
		inla_error_general(msg);
	}

	return value;
}
double inla_eval_table(char *expression, double *xval)
{
	GMRFLib_spline_tp *s;
	double *table = NULL, *x = NULL, *y = NULL, value;
	int nelm = 0, n, i, k;

	inla_sread_doubles_q(&table, &nelm, expression);
	assert(GSL_IS_EVEN(nelm));
	assert(nelm >= 4);

	n = nelm / 2;
	x = Calloc(n, double);
	y = Calloc(n, double);

	for (i = k = 0; i < n; i++, k += 2) {
		x[i] = table[k + 0];
		y[i] = table[k + 1];

		if (debug) {
			printf("table %d: %g %g\n", i, x[i], y[i]);
		}
	}

	s = inla_spline_create(x, y, n);
	value = inla_spline_eval(*xval, s);

	if (ISNAN(value)) {
		char *msg;
		GMRFLib_sprintf(&msg, "table-prior returns NAN. Argument is %g but prior is defined on [%g,%g] only.",
				*xval, s->xmin, s->xmax);
		inla_error_general(msg);
		exit(1);
	}
		
	if (0) {
		static int first = 1;
		if (first){
			FILE *fp = fopen("compare-priors.txt", "w");
			double xx;
			for(xx = -9;  xx < 9; xx += 0.01){
				double p1, p2;
				p1 = inla_spline_eval(xx, s);
				p2 = priorfunc_jeffreys_df_student_t(&xx, NULL);
				fprintf(fp, "%g %g %g\n", xx, p1, p2);
			}
			fclose(fp);
		}
		first = 0;
	}

	inla_spline_free(s);
	Free(x);
	Free(y);
	Free(table);

	return value;
}

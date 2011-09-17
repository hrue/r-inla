
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
double Gamma(double arg);
double LogGamma(double arg);
double Return(double v);
double Not(double v);
void OnError(muParserHandle_t hParser);
int set_free_variable(muParserHandle_t a_hParser, double *x);
muFloat_t *AddVariable(const muChar_t * a_szName, void *pUserData);

double Gamma(double arg)
{
	return exp(gsl_sf_lngamma(arg));
}

double LogGamma(double arg)
{
	return gsl_sf_lngamma(arg);
}

double Return(double v)
{
	return v;
}

double Not(double v)
{
	return v == 0.0;
}

void OnError(muParserHandle_t hParser)
{
	fprintf(stderr, "\n\n\nEval: Error while parsing expression:\n");
	fprintf(stderr, "----------\n");
	fprintf(stderr, "Message   :  \"%s\"\n", mupGetErrorMsg(hParser));
	fprintf(stderr, "Token     :    \"%s\"\n", mupGetErrorToken(hParser));
	fprintf(stderr, "Position  : %1d\n", mupGetErrorPos(hParser));
	fprintf(stderr, "Error Code:     %1d\n", mupGetErrorCode(hParser));
	exit(1);
}

muFloat_t *AddVariable(const muChar_t * a_szName, void *pUserData)
{
	eval_keep_vars_tp **aa = (eval_keep_vars_tp **) pUserData;
	if (*aa == NULL) {
		*aa = Calloc(1, eval_keep_vars_tp);
	}
	eval_keep_vars_tp *a = *aa;

	if (a->n >= a->n_alloc){
		a->n_alloc += 16;
		if (a->n_alloc == 0){
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
		mupSetErrorHandler(hParser, OnError);

		mupSetArgSep(hParser, ';');
		mupSetDecSep(hParser, '.');
		mupSetThousandsSep(hParser, 0);
		mupDefineConst(hParser, "pi", M_PI);
		mupDefineInfixOprt(hParser, "!", Not, 0);
		mupDefineFun1(hParser, "return", Return, 1);
		mupDefineFun1(hParser, "gamma", Gamma, 1);
		mupDefineFun1(hParser, "lgamma", LogGamma, 1);
		mupDefineFun1(hParser, "log", log, 1);
		mupDefineFun1(hParser, "log10", log10, 1);
		mupDefineFun1(hParser, "ln", log, 1);
		mupDefineFun2(hParser, "pow", pow, 1);

		eval_keep_vars_tp *keep_vars;

		keep_vars = Calloc(1, eval_keep_vars_tp);
		keep_vars->default_value = *x;
		mupSetVarFactory(hParser, AddVariable, (void *) &keep_vars);
		mupSetExpr(hParser, (muChar_t *)expression);
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

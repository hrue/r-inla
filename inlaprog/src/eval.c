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

#include <muParserDLL.h>

#include "GMRFLib/GMRFLib.h"
#include "eval.h"

/* 
   This is the interface to the muparser-library (http://muparser.sourceforge.net).  This code is buildt upon example2.c in the muparser/samples directory.
 */
    
typedef struct
{
	muFloat_t **value;
	char **name;
	int n;
}
	eval_keep_vars_tp;

/* 
   local functions...
 */
double Gamma(double arg);
double LogGamma(double arg);
double Log(double a_fVal);
double Return(double v);
double Not(double v);
void OnError(muParserHandle_t hParser);
char *get_free_variable(muParserHandle_t a_hParser);
muFloat_t* AddVariable(const muChar_t* a_szName, void *pUserData);

double Gamma(double arg) 
{
	return exp(lgamma(arg));
}

double LogGamma(double arg) 
{
	return lgamma(arg);
}

double Log(double a_fVal) 
{
	return log(a_fVal);
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
	fprintf(stderr, "\n\n\nError while parsing expression:\n");
	fprintf(stderr, "----------\n");
	fprintf(stderr, "Message   :  \"%s\"\n", mupGetErrorMsg(hParser));
	fprintf(stderr, "Token     :    \"%s\"\n", mupGetErrorToken(hParser));
	fprintf(stderr, "Position  : %1d\n", mupGetErrorPos(hParser));
	fprintf(stderr, "Error Code:     %1d\n", mupGetErrorCode(hParser));
	exit(1);
}  

char *get_free_variable(muParserHandle_t a_hParser) 
{
	muInt_t iNumVar = mupGetExprVarNum(a_hParser), i = 0;

	if (iNumVar == 0)
		return NULL;

	int num_nans = 0;
	char *name = NULL;
	for (i = 0; i < iNumVar; ++i)
	{
		const muChar_t *szName = 0;
		muFloat_t * pVar = 0;
		mupGetExprVar(a_hParser, i, &szName, &pVar);
		if (isnan(*pVar)){
			name = strdup(szName);
			num_nans++;
		}
	}
	if (num_nans != 1){
		fprintf(stderr, "\n\nError: There are %1d free variables in the expression. Only one is allowed.\n\n", num_nans);
		exit(1);
	}

	return name;
}  

muFloat_t* AddVariable(const muChar_t* a_szName, void *pUserData)
{
	eval_keep_vars_tp **aa = (eval_keep_vars_tp **) pUserData;
	if (*aa == NULL) {
		*aa = Calloc(1, eval_keep_vars_tp);
	}
	eval_keep_vars_tp *a = *aa;

	a->value = Realloc(a->value, a->n+1,  muFloat_t *);
	a->name = Realloc(a->name, a->n+1,  char *);
	a->value[a->n] = Calloc(1, muFloat_t);
	a->value[a->n][0] = NAN;
	a->name[a->n] = strdup(a_szName);
	a->n++;

	return &(a->value[a->n-1][0]);
}

double inla_eval(char *expression, double x) 
{
	muParserHandle_t hParser;

	hParser = mupCreate();	
	mupSetErrorHandler(hParser, OnError);
	 
	mupSetArgSep(hParser, ';');
	mupSetDecSep(hParser, '.');
	mupSetThousandsSep(hParser, 0);
	mupDefineConst(hParser, "pi", M_PI);
	mupDefineInfixOprt(hParser, "!", Not, 0);
	mupDefineFun1(hParser, "return", Return, 0);
	mupDefineFun1(hParser, "log", Log, 1);
	mupDefineFun1(hParser, "gamma", Gamma, 1);
	mupDefineFun1(hParser, "lgamma", LogGamma, 1);

	eval_keep_vars_tp *keep_vars = NULL;
	mupSetVarFactory(hParser, AddVariable, (void *) &keep_vars);

	mupSetExpr(hParser, (muChar_t *)expression);

	char *name_x = get_free_variable(hParser);
	mupDefineVar(hParser, name_x, &x);
	muFloat_t value = mupEval(hParser);

	mupRelease(hParser);
	Free(name_x);
	if (keep_vars) {
		int i;
		for(i=0; i<keep_vars->n; i++){

			printf("Free %s = %g\n", keep_vars->name[i], (double) keep_vars->value[i][0]);
			Free(keep_vars->value[i]);
			Free(keep_vars->name[i]);
		}
		Free(keep_vars->value);
		Free(keep_vars->name);
		Free(keep_vars);
	}
	
	return (double) value;
}

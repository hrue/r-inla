
/* example-wa.c
 * 
 * Copyright (C) 2001 Havard Rue
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
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA.
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

#include <string.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif

static const char RCSId[] = "$Id: example-wa.c,v 1.12 2008/10/29 16:32:25 hrue Exp $";

#include "GMRFLib/GMRFLib.h"

#define DATAFILE "wa-example.dat"
#define N 50
#define NSAMPLES 10

double wafunc(int i, int j, void *arg);

double wafunc(int i, int j, void *arg)
{
	/*
	 * return the weights w_ij, including the precision sqrt(h) 
	 */
	double val, h;
	void **args;
	int order;

	args = (void **) arg;
	order = *(int *) args[0];
	h = *(double *) args[1];

	if (j > i || i - j < order) {
		return 0.0;
	}

	switch (order) {
	case 1:
		return ((j == i || j == i - 1) ? sqrt(h) : 0.0);
	case 2:
		switch (i - j) {
		case 0:
			val = 1.0;
			break;
		case 1:
			val = 2.0;
			break;
		case 2:
			val = -1.0;
			break;
		default:
			val = 0.0;
		}
		return val * sqrt(h);
	default:
		abort();
	}
	return 0.0;
}
int main(int argc, char **argv)
{
	int seed, i, k, order;
	double h = 100.0, *x, *c, *b, *y, *samples, noise_variance;
	FILE *fp;
	void *args[2];
	GMRFLib_graph_tp *wagraph;
	GMRFLib_wa_problem_tp *waproblem;
	GMRFLib_problem_tp *problem;

	if (argc < 4) {
		printf("Usage: %s SEED ORDER H\n", *argv);
		return 0;
	}

	/*
	 * get the seed and init the random generator in GMRFLib, read other arguments 
	 */
	seed = atoi(argv[1]);
	GMRFLib_uniform_init((unsigned long int) seed);

	order = atoi(argv[2]);
	h = atof(argv[3]);

	/*
	 * setup arguments to wafunc 
	 */
	args[0] = (void *) &order;
	args[1] = (void *) &h;

	/*
	 * create wa-problem, note that the graph, Qfunc and Qfunc_arg is then provided in 'waproblem'! 
	 */
	GMRFLib_make_linear_graph(&wagraph, N, order, 0);
	GMRFLib_init_wa_problem(&waproblem, wagraph, wafunc, (void *) args);

	/*
	 * setup space 
	 */
	x = Calloc(N, double);
	c = Calloc(N, double);
	b = Calloc(N, double);
	y = Calloc(N, double);
	samples = Calloc(N * NSAMPLES, double);

	/*
	 * read data 
	 */
	noise_variance = SQR(0.25);
	fp = fopen(DATAFILE, "r");
	for (i = 0; i < N; i++)
		fscanf(fp, "%lf", &y[i]);
	fclose(fp);

	/*
	 * setup the posterior for x|y 
	 */
	for (i = 0; i < N; i++) {
		b[i] = y[i] / noise_variance;
		c[i] = 1.0 / noise_variance;
	}

	/*
	 * compute the conditonal mean and NSAMPLES samples 
	 */

	GMRFLib_init_problem(&problem, x, b, c, NULL, waproblem->graph, waproblem->Qfunc, waproblem->Qfunc_arg, NULL, NULL, GMRFLib_NEW_PROBLEM);
	for (k = 0; k < NSAMPLES; k++) {
		GMRFLib_sample(problem);
		memcpy(&samples[k * N], problem->sample, N * sizeof(double));
	}

	for (i = 0; i < N; i++) {
		printf("%d %f %f", i, y[i], problem->sub_mean[i]);
		for (k = 0; k < NSAMPLES; k++) {
			printf(" %f", samples[k * N + i]);
		}
		printf("\n");
	}

	GMRFLib_free_problem(problem);
	GMRFLib_free_wa_problem(waproblem);
	GMRFLib_free_graph(wagraph);
	Free(x);
	Free(c);
	Free(b);
	Free(y);
	Free(samples);

	return 0;
}


/* example-rw.c
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

#include <assert.h>
#include <string.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif

static const char RCSId[] = "$Id: example-rw.c,v 1.12 2008/10/29 16:32:25 hrue Exp $";

#include "GMRFLib/GMRFLib.h"

#define DATAFILE "wa-example.dat"
#define N 50
#define NSAMPLES 10

int main(int argc, char **argv)
{
	/*
	 * solve the same problem as in example-wa.c now using the routines in rw.c 
	 */

	int seed, i, k, order;
	double h = 100.0, *x = NULL, *c = NULL, *b = NULL, *y = NULL, *samples, noise_variance;
	FILE *fp;
	GMRFLib_graph_tp *graph;
	GMRFLib_problem_tp *problem;
	GMRFLib_rwdef_tp *rwdef = NULL;
	GMRFLib_crwdef_tp *crwdef = NULL;

	GMRFLib_smtp = GMRFLib_SMTP_BAND;		       /* use band alg here */

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

	assert(order == 1 || order == 2);
	assert(h > 0.0);

	if (order == 1) {
		/*
		 * use RW1 
		 */
		rwdef = Calloc(1, GMRFLib_rwdef_tp);	       /* also ensure that all is NULL */
		rwdef->n = N;
		rwdef->order = 1;
		rwdef->cyclic = 0;
		rwdef->prec = &h;

		GMRFLib_make_rw_graph(&graph, rwdef);
	} else if (order == 2) {
		/*
		 * use CRW2 
		 */
		crwdef = Calloc(1, GMRFLib_crwdef_tp);	       /* also ensure that all is NULL */
		crwdef->n = N;
		crwdef->order = 2;
		crwdef->prec = &h;
		crwdef->layout = GMRFLib_CRW_LAYOUT_BLOCK;
		GMRFLib_make_crw_graph(&graph, crwdef);
	}

	/*
	 * setup space 
	 */
	if (order == 1) {
		x = Calloc(N, double);
		c = Calloc(N, double);
		b = Calloc(N, double);
		y = Calloc(N, double);
	} else if (order == 2) {
		x = Calloc(2 * N, double);
		c = Calloc(2 * N, double);
		b = Calloc(2 * N, double);
		y = Calloc(2 * N, double);
	}

	samples = Calloc(N * NSAMPLES, double);

	/*
	 * read data 
	 */
	noise_variance = SQR(0.2500);
	fp = fopen(DATAFILE, "r");
	for (i = 0; i < N; i++) {
		fscanf(fp, "%lf", &y[i]);
	}
	fclose(fp);

	/*
	 * setup the posterior for x|y. this is ok also for CRW2 since we use a BLOCKED layout. 
	 */
	for (i = 0; i < N; i++) {
		b[i] = y[i] / noise_variance;
		c[i] = 1.0 / noise_variance;
	}

	/*
	 * compute the conditonal mean and NSAMPLES samples 
	 */

	GMRFLib_init_problem(&problem, x, b, c, NULL, graph,
			     (order == 1 ? GMRFLib_rw : GMRFLib_crw), (order == 1 ? (void *) rwdef : (void *) crwdef), NULL, NULL,
			     GMRFLib_NEW_PROBLEM);
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

	return 0;
}

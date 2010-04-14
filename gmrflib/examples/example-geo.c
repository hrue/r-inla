
/* example-geo.c
 * 
 * Copyright (C) 2005 Havard Rue
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
#include <stddef.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif

static const char RCSId[] = "$Id: example-geo.c,v 1.22 2010/03/12 12:24:40 hrue Exp $";

#include "GMRFLib/GMRFLib.h"

int main(int argc, char **argv)
{
	GMRFLib_geo_problem_tp *geo_problem;
	int name, neigh, nrow, ncol, cyclic_flag;
	double param, range, prec, *marg_var;
	int i, j, idx;
	GMRFLib_problem_tp *problem;
	FILE *fp;

	if (argc < 2) {
		fprintf(stderr, "Usage: %s SEED\n", argv[0]);
		exit(1);
	}

	if ((fp = fopen("available-geo-coefs.dat", "w"))) {
		printf("print available geo-coefs to file\n");
		GMRFLib_print_geo_coefs(fp);
		fclose(fp);
	}

	if (0) {
		/*
		 * check available cooofs 
		 */

		if (GMRFLib_is_geo_coefs(GMRFLib_CORTP_MATERN, 2, 1.0, 10.0) == GMRFLib_SUCCESS) {
			printf("(GMRFLib_CORTP_MATERN, 2, 1.0, 10.0) exists: OK\n");
		} else {
			printf("(GMRFLib_CORTP_MATERN, 2, 1.0, 10.0) does not exists: FATAL ERROR\n");
			abort();
		}

		if (GMRFLib_is_geo_coefs(GMRFLib_CORTP_MATERN, 2, 1.11, 11.11) == GMRFLib_SUCCESS) {
			printf("(GMRFLib_CORTP_MATERN, 2, 1.11, 11.11) exists: weird?\n");
		} else {
			printf("(GMRFLib_CORTP_MATERN, 2, 1.11, 11.11) does not exists: OK\n");
		}
	}

	name = GMRFLib_CORTP_MATERN;
	neigh = 2;
	prec = 1.0;
	param = 1.5;
	range = 16.0;
	nrow = 10;
	ncol = 10;
	cyclic_flag = 0;

	GMRFLib_init_geo_problem(&geo_problem, name, neigh, param, range, nrow, ncol, &prec, cyclic_flag);

	/*
	 * print the Q-matrix 
	 */
	if (1) {
		fp = fopen("Q.dat", "w");
		printf("write [Q.dat]\n");
		for (i = 0; i < geo_problem->graph->n; i++) {
			fprintf(fp, "%d %d %.16f\n", i, i, geo_problem->Qfunc(i, i, geo_problem->Qfunc_arg));
			for (j = 0; j < geo_problem->graph->nnbs[i]; j++) {
				fprintf(fp, "%d %d %.16f\n", i, geo_problem->graph->nbs[i][j],
					geo_problem->Qfunc(i, geo_problem->graph->nbs[i][j], geo_problem->Qfunc_arg));
			}
		}
		fclose(fp);
	}

	/*
	 * sample a geo-problem 
	 */

	if (1) {
		GMRFLib_uniform_init((unsigned long int) atol(argv[1]));	// init the GMRFLib's random generator

		GMRFLib_init_problem(&problem, NULL, NULL, NULL, NULL, geo_problem->graph,
				     geo_problem->Qfunc, geo_problem->Qfunc_arg, NULL, NULL, GMRFLib_NEW_PROBLEM);

		GMRFLib_sample(problem);
		fp = fopen("sample.dat", "w");
		printf("write [sample.dat]\n");
		for (i = 0; i < nrow; i++) {
			for (j = 0; j < ncol; j++) {
				GMRFLib_lattice2node(&idx, i, j, nrow, ncol);
				fprintf(fp, " %lf", problem->sample[idx]);
			}
			fprintf(fp, "\n");
		}
		fclose(fp);

		GMRFLib_Qinv(problem, GMRFLib_QINV_ALL);
		fp = fopen("variances.dat", "w");
		printf("write [variances.dat]\n");
		for (i = 0; i < nrow; i++) {
			for (j = 0; j < ncol; j++) {
				GMRFLib_lattice2node(&idx, i, j, nrow, ncol);
				marg_var = GMRFLib_Qinv_get(problem, idx, idx);
				fprintf(fp, " %lf", *marg_var);
			}
			fprintf(fp, "\n");
		}
		fclose(fp);
	}

	if (1) {
		/*
		 * illustrate the revision of an GMRFLib_geo_problem_tp-object 
		 */

		name = GMRFLib_CORTP_EXP;		       /* new type */
		prec = 1.0;
		param = -1.0;				       /* not in use for EXP, set to any value */
		range = 40.0;				       /* new range */

		GMRFLib_revise_geo_problem(geo_problem, name, param, range, &prec);

		/*
		 * init again reusing the graph. [not that big saving, but....] 
		 */
		GMRFLib_init_problem(&problem, NULL, NULL, NULL, NULL, geo_problem->graph,
				     geo_problem->Qfunc, geo_problem->Qfunc_arg, NULL, NULL, GMRFLib_KEEP_graph);

		GMRFLib_sample(problem);

		fp = fopen("sample-revised.dat", "w");
		printf("write [sample-revised.dat]\n");
		for (i = 0; i < nrow; i++) {
			for (j = 0; j < ncol; j++) {
				GMRFLib_lattice2node(&idx, i, j, nrow, ncol);
				fprintf(fp, " %lf", problem->sample[idx]);
			}
			fprintf(fp, "\n");
		}
		fclose(fp);
	}

	GMRFLib_free_geo_problem(geo_problem);

	return 0;
}

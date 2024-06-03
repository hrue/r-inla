
/* param-constr.c
 * 
 * Copyright (C) 2024-2024 Havard Rue
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

#define debug_(msg_)							\
	if (debug) {							\
		printf("\t[%1d] %s (%s:%1d): %s\n", omp_get_thread_num(), GMRFLib_function_name_strip(__GMRFLib_FuncName), __FILE__, __LINE__, msg_); \
	}

#define debug2_(msg_, val_)						\
	if (debug) {							\
		printf("\t[%1d] %s (%s:%1d): %s %g\n", omp_get_thread_num(), GMRFLib_function_name_strip(__GMRFLib_FuncName), __FILE__, __LINE__, msg_, val_); \
	}

int inla_parse_param_constraints(inla_tp *mb)
{
	const int debug = 1;
	for (int k = 0; k < mb->nds; k++) {
		Data_section_tp *ds = &(mb->data_sections[k]);
		if (!(ds->data_observations.sem_B)) {
			continue;
		}

		int dim = ds->data_observations.sem_dim;
		if (debug) {
			printf("parse section %d dim %d\n", k, dim);
		}

		ds->data_observations.sem_B_ptr = Calloc(ISQR(dim), double **);
		ds->data_observations.sem_B_map = Calloc(ISQR(dim), map_func_tp *);
		ds->data_observations.sem_B_map_arg = Calloc(ISQR(dim), void *);
		char **B = ds->data_observations.sem_B;	       /* shorthand */

		for (int i = 0; i < ISQR(dim); i++) {
			if (debug) {
				printf("search for name B[%1d] = [%s]\n", i, B[i]);
			}

			if (strcasecmp(B[i], SEM_NO_VAR) != 0) {
				int j = find_tag(mb, B[i]);
				if (j < 0) {
					char *err = NULL;
					GMRFLib_sprintf(&err, "SEM: search for model component B[%1d] = [%s]: not found\n", i, B[i]);
					inla_error_general(err);
					exit(1);
				}
				if (mb->f_id[j] != F_COPY) {
					char *err = NULL;
					GMRFLib_sprintf(&err, "SEM: model component [%s] not COPY\n", mb->f_tag[j]);
					inla_error_general(err);
					exit(1);
				}
				ds->data_observations.sem_B_ptr[i] = mb->f_theta[j][0];
				ds->data_observations.sem_B_map[i] = mb->f_theta_map[j][0];
				ds->data_observations.sem_B_map_arg[i] = mb->f_theta_map_arg[j][0];
			} else {
				double **one = NULL;
				HYPER_NEW(one, 1.0);
				ds->data_observations.sem_B_ptr[i] = one;
				ds->data_observations.sem_B_map[i] = map_one;
				ds->data_observations.sem_B_map_arg[i] = NULL;
			}

			if (debug) {
				double bvalue = ds->data_observations.sem_B_map[i] (ds->data_observations.sem_B_ptr[i][0][0],
										    MAP_FORWARD,
										    ds->data_observations.sem_B_map_arg[i]);
				printf("B[%1d] :  A = %g  B = %g\n", i, ds->data_observations.sem_A[i], bvalue);
			}
		}
	}

	return 0;
}

double inla_eval_param_constraint(int thread_id, Data_section_tp *ds)
{
	typedef struct {
		double *theta;
		double value;
		size_t dim;
		int idx;
	} cache_tp;

	const int debug = 0;
	int idx = ds->data_observations.sem_idx;
	size_t dim = ds->data_observations.sem_dim;
	int in_cache = 1;

	if (!(ds->data_observations.sem_cache)) {
#pragma omp critical  (Name_63a8d38d09184295d229c3a47bcbba98187ddc96)
		if (!(ds->data_observations.sem_cache)) {
			cache_tp **c = Calloc(GMRFLib_CACHE_LEN(), cache_tp *);
			ds->data_observations.sem_cache = (void *) c;
		}
	}

	int cache_idx = -1;
	GMRFLib_CACHE_SET_ID(cache_idx);
	cache_tp *c = ((cache_tp **) ds->data_observations.sem_cache)[cache_idx];

	if (!c) {
#pragma omp critical (Name_518e2a01a0639bcfa1a50ef2e26dd6d69506dacf)
		if (!c) {
			debug_("INIT CACHE");
			cache_tp *cc = Calloc(1, cache_tp);
			cc->dim = dim;
			cc->theta = Calloc(ISQR(dim), double);
			in_cache = 0;
			c = ((cache_tp **) ds->data_observations.sem_cache)[cache_idx] = cc;
		}
	}

	if (c->dim != dim) {
#pragma omp critical (Name_3f524981b6edda449347e31549f5e8b399411d79)
		if (c->dim != dim) {
			debug_("DIM CHANGE IN CACHE");
			if (dim > c->dim) {
				c->theta = Realloc(c->theta, IMAX(ISQR(dim), ISQR(c->dim)), double);
			}
			in_cache = 0;
			c->dim = dim;
		}
	}

	if (c->idx != idx) {
#pragma omp critical (Name_bd076b6eeef2e2fce722490f1be7ba6118874702)
		if (c->idx != idx) {
			debug_("IDX CHANGE IN CACHE");
			in_cache = 0;
			c->idx = idx;
		}
	}

	if (in_cache) {
		for (int k = 0; k < ISQR(dim); k++) {
			if (c->theta[k] != ds->data_observations.sem_B_ptr[k][thread_id][0]) {
				in_cache = 0;
				break;
			}
		}
		if (in_cache) {
			return (c->value);
		}
	}

	gsl_matrix *B = gsl_matrix_calloc(dim, dim);
	gsl_matrix *S = gsl_matrix_calloc(dim, dim);
	for (size_t j = 0, k = 0; j < dim; j++) {
		for (size_t i = 0; i < dim; i++) {
			double Bval, Aval, ABval;
			c->theta[k] = ds->data_observations.sem_B_ptr[k][thread_id][0];
			Bval = ds->data_observations.sem_B_map[k] (c->theta[k], MAP_FORWARD, ds->data_observations.sem_B_map_arg[k]);
			Aval = ds->data_observations.sem_A[k];
			ABval = Aval * Bval;
			gsl_matrix_set(B, i, j, (i == j ? 1.0 - ABval : -ABval));
			k++;
		}
	}

	if (0 && debug) {
		debug_("I-B\n");
		GMRFLib_printf_gsl_matrix(stdout, B, " %g");
	}


	GMRFLib_gsl_mmt(B, B, S);
	GMRFLib_gsl_spd_inverse(S);			       /* S <- solve(B %*% t(B)) */
	c->value = 1.0 / gsl_matrix_get(S, idx, idx);	       /* return the marginal precision */

	debug2_("COMPUTE NEW VALUE FOR PRECISION IN CACHE", c->value);

	gsl_matrix_free(B);
	gsl_matrix_free(S);

	return c->value;
}
#undef debug_
#undef debug2_

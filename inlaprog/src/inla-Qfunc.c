
/* inla-Qfunc.c
 * 
 * Copyright (C) 2007-2024 Havard Rue
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 *  * This program is distributed in the hope that it will be useful, but
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

double Qfunc_bym(int thread_id, int i, int j, double *UNUSED(values), void *arg)
{
	if (j < 0) {
		return NAN;
	}

	/*
	 * the first half is the `linear predictor', the second half is the spatial term. x = [eta, z], eta|.. ~ N(z, prec*I), z ~ GMRF(prec) 
	 */

	inla_bym_Qfunc_arg_tp *a = (inla_bym_Qfunc_arg_tp *) arg;
	int n = a->n;
	double prec_iid = map_precision(a->log_prec_iid[thread_id][0], MAP_FORWARD, NULL);

	if (IMAX(i, j) < n) {
		/*
		 * the iid term 
		 */
		return prec_iid;
	}
	if (IMIN(i, j) >= n) {
		/*
		 * the spatial term + I*prec_iid 
		 */
		return (i == j ? prec_iid : 0.0) + Qfunc_besag(thread_id, i - n, j - n, NULL, a->besag_arg);
	} else {
		/*
		 * the cross-term which is -prec_iid 
		 */
		return -prec_iid;
	}

	/*
	 * should not happen 
	 */
	GMRFLib_ASSERT(0 == 1, GMRFLib_ESNH);
	return 0.0;
}

double Qfunc_bym2(int thread_id, int i, int j, double *UNUSED(values), void *arg)
{
	if (j < 0) {
		return NAN;
	}

	inla_bym2_Qfunc_arg_tp *a = (inla_bym2_Qfunc_arg_tp *) arg;
	int n = a->n;
	double prec = map_precision(a->log_prec[thread_id][0], MAP_FORWARD, NULL);
	double phi = map_probability(a->logit_phi[thread_id][0], MAP_FORWARD, NULL);
	if (IMAX(i, j) < n) {
		return prec / (1.0 - phi);
	}
	if (IMIN(i, j) >= n) {
		return (i == j ? phi / (1.0 - phi) : 0.0) + Qfunc_besag(thread_id, i - n, j - n, NULL, a->besag_arg);
	}
	return -sqrt(phi * prec) / (1.0 - phi);
}

double Qfunc_rw2diid(int thread_id, int i, int j, double *UNUSED(values), void *arg)
{
	if (j < 0) {
		return NAN;
	}

	inla_rw2diid_Qfunc_arg_tp *a = (inla_rw2diid_Qfunc_arg_tp *) arg;
	int n = a->n;
	double prec = map_precision(a->log_prec[thread_id][0], MAP_FORWARD, NULL);
	double phi = map_probability(a->logit_phi[thread_id][0], MAP_FORWARD, NULL);
	if (IMAX(i, j) < n) {
		return prec / (1.0 - phi);
	}
	if (IMIN(i, j) >= n) {
		return (i == j ? phi / (1.0 - phi) : 0.0) + GMRFLib_rw2d(thread_id, i - n, j - n, NULL, a->rw2ddef);
	}
	return -sqrt(phi * prec) / (1.0 - phi);
}

double Qfunc_group(int thread_id, int i, int j, double *UNUSED(values), void *arg)
{
	// spell out map_precision(...)
	if (j < 0) {
		return NAN;
	}

	inla_group_def_tp *a = (inla_group_def_tp *) arg;
	double val = 0.0, fac = 0.0, rho = 0.0, prec = 0.0;

	int n = a->N;					       /* this is the size before group */
	int ngroup = a->ngroup;

	div_t ii = div(i, n);
	div_t jj = div(j, n);
	
	int igroup = ii.quot;
	int jgroup = jj.quot;
	int irem = ii.rem;
	int jrem = jj.rem;
	
	switch (a->type) {
	case G_EXCHANGEABLE: 
	{
		rho = map_group_rho(a->group_rho_intern[thread_id][0], MAP_FORWARD, (void *) &(a->ngroup));
		if (igroup == jgroup) {
			fac = -((ngroup - 2.0) * rho + 1.0) / ((rho - 1.0) * ((ngroup - 1.0) * rho + 1.0));
		} else {
			fac = rho / ((rho - 1.0) * ((ngroup - 1.0) * rho + 1.0));
		}
	}
	break;

	case G_EXCHANGEABLE_POS:
	{
		rho = map_probability(a->group_rho_intern[thread_id][0], MAP_FORWARD, (void *) &(a->ngroup));
		if (igroup == jgroup) {
			fac = -((ngroup - 2.0) * rho + 1.0) / ((rho - 1.0) * ((ngroup - 1.0) * rho + 1.0));
		} else {
			fac = rho / ((rho - 1.0) * ((ngroup - 1.0) * rho + 1.0));
		}
	}
	break;

	case G_AR1:
	{
		rho = map_rho(a->group_rho_intern[thread_id][0], MAP_FORWARD, NULL);
		if (igroup == jgroup) {
			if (!(a->cyclic) && (igroup == 0 || igroup == ngroup - 1)) {
				fac = 1.0 / (1.0 - SQR(rho));
			} else {
				fac = (1.0 + SQR(rho)) / (1.0 - SQR(rho));
			}
		} else {
			fac = -rho / (1.0 - SQR(rho));
		}
	}
	break;

	case G_AR:
	{
		fac = Qfunc_ar(thread_id, igroup, jgroup, NULL, (void *) a->ardef);
	}
	break;

	case G_RW1:
	case G_RW2:
	{
		//prec = map_precision(a->group_prec_intern[thread_id][0], MAP_FORWARD, NULL);
		prec = exp(a->group_prec_intern[thread_id][0]);
		if (a->crwdef) {
			fac = prec * GMRFLib_crw(thread_id, igroup, jgroup, NULL, (void *) (a->crwdef));
		} else {
			fac = prec * GMRFLib_rw(thread_id, igroup, jgroup, NULL, (void *) (a->rwdef));
		}
	}
	break;

	case G_BESAG:
	{
		//prec = map_precision(a->group_prec_intern[thread_id][0], MAP_FORWARD, NULL);
		prec = exp(a->group_prec_intern[thread_id][0]);
		fac = prec * Qfunc_besag(thread_id, igroup, jgroup, NULL, (void *) (a->besagdef));
	}
	break;

	case G_IID:
	{
		if (igroup == jgroup) {
			//fac = map_precision(a->group_prec_intern[thread_id][0], MAP_FORWARD, NULL);
			fac = exp(a->group_prec_intern[thread_id][0]);
		} else {
			fac = 0.0;
		}
	}
	break;

	default:
		assert(0 == 1);
	}

	val = a->Qfunc(thread_id, irem, jrem, NULL, a->Qfunc_arg) * fac;
	return val;
}

double Qfunc_generic1(int thread_id, int i, int j, double *UNUSED(values), void *arg)
{
	if (j < 0) {
		return NAN;
	}

	inla_generic1_tp *a = (inla_generic1_tp *) arg;
	double beta = map_probability(a->beta[thread_id][0], MAP_FORWARD, NULL);
	double prec = map_precision(a->log_prec[thread_id][0], MAP_FORWARD, NULL);

	return prec * ((i == j ? 1.0 : 0.0) - (beta / a->max_eigenvalue) * a->tab->Qfunc(thread_id, i, j, NULL, a->tab->Qfunc_arg));
}

double Qfunc_generic2(int thread_id, int i, int j, double *UNUSED(values), void *arg)
{
	if (j < 0) {
		return NAN;
	}

	/*
	 * [u,v], where u = v + noise. h2 = 1/prec_v / (1/prec_v + 1/prev_u). h2 in (0,1). 
	 */

	inla_generic2_tp *a = (inla_generic2_tp *) arg;

	double prec_cmatrix = map_precision(a->log_prec[thread_id][0], MAP_FORWARD, NULL),
	    h2 = map_probability(a->h2_intern[thread_id][0], MAP_FORWARD, NULL), prec_unstruct;
	int n = a->n;

	prec_unstruct = h2 / (1.0 - h2) * prec_cmatrix;

	if (i == j) {

		if (i < n) {
			/*
			 * the sum 
			 */
			return prec_unstruct;
		} else {
			/*
			 * cmatrix 
			 */
			return prec_unstruct + prec_cmatrix * a->tab->Qfunc(thread_id, i - n, i - n, NULL, a->tab->Qfunc_arg);
		}
	} else {
		if (IMIN(i, j) >= n) {
			/*
			 * inside the Cmatrix 
			 */
			return prec_cmatrix * a->tab->Qfunc(thread_id, i - n, j - n, NULL, a->tab->Qfunc_arg);
		} else {
			/*
			 * cross term between cmatrix and the sum 
			 */
			return -prec_unstruct;
		}
	}
	GMRFLib_ASSERT(0 == 1, GMRFLib_ESNH);
	return 0.0;
}

double Qfunc_generic3(int thread_id, int i, int j, double *UNUSED(values), void *arg)
{
	if (j < 0) {
		return NAN;
	}

	inla_generic3_tp *a = (inla_generic3_tp *) arg;
	double prec_common, prec, val = 0.0;
	int k, same = (i == j);

	for (k = 0; k < a->m; k++) {
		if (same || GMRFLib_graph_is_nb(i, j, a->g[k])) {
			prec = map_precision(a->log_prec[k][thread_id][0], MAP_FORWARD, NULL);
			val += prec * a->tab[k]->Qfunc(thread_id, i, j, NULL, a->tab[k]->Qfunc_arg);
		}
	}
	prec_common = map_precision(a->log_prec[GENERIC3_MAXTHETA - 1][thread_id][0], MAP_FORWARD, NULL);
	val *= prec_common;

	return (val);
}

double Qfunc_replicate(int thread_id, int i, int j, double *UNUSED(values), void *arg)
{
	if (j < 0) {
		return NAN;
	}

	int ii, jj;
	inla_replicate_tp *a = (inla_replicate_tp *) arg;

	ii = MOD(i, a->n);
	jj = MOD(j, a->n);

	return a->Qfunc(thread_id, ii, jj, NULL, a->Qfunc_arg);
}

double Qfunc_z(int thread_id, int i, int j, double *UNUSED(values), void *arg)
{
	if (j < 0) {
		return NAN;
	}

	inla_z_arg_tp *a = (inla_z_arg_tp *) arg;
	double value = 0.0;

	if (i == j || GMRFLib_graph_is_nb(i, j, a->graph_A)) {
		value += a->Qfunc_A->Qfunc(thread_id, i, j, NULL, a->Qfunc_A->Qfunc_arg);
	}
	if (i == j || GMRFLib_graph_is_nb(i, j, a->graph_B)) {
		/*
		 * doit like this, as most of the elements in B are zero
		 */
		double q = a->Qfunc_B->Qfunc(thread_id, i, j, NULL, a->Qfunc_B->Qfunc_arg);
		if (q) {
			value += q * map_precision(a->log_prec[thread_id][0], MAP_FORWARD, NULL);
		}
	}
	return value;
}

double Qfunc_slm(int thread_id, int i, int j, double *UNUSED(values), void *arg)
{
	if (j < 0) {
		return NAN;
	}

	inla_slm_arg_tp *a = (inla_slm_arg_tp *) arg;
	double value = 0.0, prec, rho, rho_std;

	prec = map_precision(a->log_prec[thread_id][0], MAP_FORWARD, NULL);
	rho_std = map_probability(a->logit_rho[thread_id][0], MAP_FORWARD, NULL);
	rho = a->rho_min + rho_std * (a->rho_max - a->rho_min);

	if (0) {
		if (i == j && i == 0) {
			printf("log_prec %.12f  logit_rho %.12f\n", a->log_prec[thread_id][0], a->logit_rho[thread_id][0]);
		}
	}

	if (i == j || GMRFLib_graph_is_nb(i, j, a->graph_A1)) {
		value += prec * a->Qfunc_A1->Qfunc(thread_id, i, j, NULL, a->Qfunc_A1->Qfunc_arg);
	}
	if (i == j || GMRFLib_graph_is_nb(i, j, a->graph_A2)) {
		value += a->Qfunc_A2->Qfunc(thread_id, i, j, NULL, a->Qfunc_A2->Qfunc_arg);
	}
	if (i == j || GMRFLib_graph_is_nb(i, j, a->graph_B)) {
		value += prec * rho * a->Qfunc_B->Qfunc(thread_id, i, j, NULL, a->Qfunc_B->Qfunc_arg);
	}
	if (i == j || GMRFLib_graph_is_nb(i, j, a->graph_C)) {
		value += prec * SQR(rho) * a->Qfunc_C->Qfunc(thread_id, i, j, NULL, a->Qfunc_C->Qfunc_arg);
	}

	return value;
}

double Qfunc_rgeneric(int thread_id, int i, int j, double *values, void *arg)
{
	inla_rgeneric_tp *a = (inla_rgeneric_tp *) arg;
	int rebuild, ii, id = 0;
	const int debug = 0;

	GMRFLib_CACHE_SET_ID(id);
	rebuild = (a->param[id] == NULL || a->Q[id] == NULL);

	if (!rebuild) {
		for (ii = 0; ii < a->ntheta && !rebuild; ii++) {
			rebuild = (a->param[id][ii] != a->theta[ii][thread_id][0]);
		}
	}

	if (rebuild) {
		int *ilist = NULL, *jlist = NULL, n, len, k = 0, n_out, jj;
		double *Qijlist = NULL, *x_out = NULL;
#pragma omp critical (Name_297cd7aba8c5dafefcb1c93779913d23a945ce9e)
		{
			rebuild = (a->param[id] == NULL || a->Q[id] == NULL);
			if (!rebuild) {
				for (ii = 0; ii < a->ntheta && !rebuild; ii++) {
					rebuild = (a->param[id][ii] != a->theta[ii][thread_id][0]);
				}
			}

			if (rebuild) {
				if (debug) {
					printf("Qfunc_rgeneric: Rebuild Q-hash for id %d thread_id %d\n", id, thread_id);
				}
				if (a->Q[id]) {
					GMRFLib_free_tabulate_Qfunc(a->Q[id]);
				}
				double *a_tmp = Calloc(a->ntheta, double);
				for (jj = 0; jj < a->ntheta; jj++) {
					a_tmp[jj] = a->theta[jj][thread_id][0];
					if (debug) {
						printf("\ttheta[%1d] %.12f\n", jj, a_tmp[jj]);
					}
				}

				if (debug) {
					printf("\tCall rgeneric\n");
				}
				inla_R_rgeneric(&n_out, &x_out, R_GENERIC_Q, a->model, &(a->ntheta), a_tmp);
				if (debug) {
					printf("\tReturn from rgeneric with n_out= %1d\n", n_out);
				}
				assert(n_out >= 2);

				assert(a->graph);
				if ((int) x_out[0] == -1) {
					// optimized output
					k = 1;
					len = (int) x_out[k++];
					assert(len == a->len_list);
					n = a->graph->n;
					GMRFLib_tabulate_Qfunc_from_list2(&(a->Q[id]), a->graph, a->len_list, a->ilist, a->jlist, &(x_out[k]), n,
									  NULL);
				} else {
					n = (int) x_out[k++];
					len = (int) x_out[k++];

					// we can overlay these arrays to avoid allocating new ones, since x_out is double
					ilist = (int *) &(x_out[k]);
					jlist = (int *) &(x_out[k + len]);
					Qijlist = (double *) &(x_out[k + 2 * len]);
					for (jj = 0; jj < len; jj++) {
						ilist[jj] = (int) x_out[k + jj];
						jlist[jj] = (int) x_out[k + len + jj];
					}

					GMRFLib_tabulate_Qfunc_from_list2(&(a->Q[id]), a->graph, len, ilist, jlist, Qijlist, n, NULL);
					assert(a->graph->n == a->n);
				}
				Free(x_out);

				if (a->param[id]) {
					Memcpy(a->param[id], a_tmp, a->ntheta * sizeof(double));
					Free(a_tmp);
				} else {
					a->param[id] = a_tmp;
				}
				if (debug) {
					printf("\tRebuild for id %1d done\n", id);
				}
			}
		}
	}

	if (j >= 0) {
		return (a->Q[id]->Qfunc(thread_id, i, j, NULL, a->Q[id]->Qfunc_arg));
	} else {
		return (a->Q[id]->Qfunc(thread_id, i, j, values, a->Q[id]->Qfunc_arg));
	}
}

double Qfunc_cgeneric(int thread_id, int i, int j, double *values, void *arg)
{
	inla_cgeneric_tp *a = (inla_cgeneric_tp *) arg;
	int rebuild, id = 0;
	const int debug = 0;

	GMRFLib_CACHE_SET_ID(id);
	rebuild = (a->param[id] == NULL || a->Q[id] == NULL);
	if (!rebuild) {
		for (int ii = 0; ii < a->ntheta && !rebuild; ii++) {
			rebuild = (a->param[id][ii] != a->theta[ii][thread_id][0]);
		}
	}

	if (rebuild) {
		int *ilist = NULL, *jlist = NULL, n, len, k = 0;
		double *Qijlist = NULL, *x_out = NULL;
		rebuild = (a->param[id] == NULL || a->Q[id] == NULL);
		if (!rebuild) {
			for (int ii = 0; ii < a->ntheta && !rebuild; ii++) {
				rebuild = (a->param[id][ii] != a->theta[ii][thread_id][0]);
				if (debug) {
					printf("\t\tid= %1d thread_id=%1d theta[%1d] %.12f\n", id, thread_id, ii, a->theta[ii][thread_id][0]);
				}
			}
		}

		if (rebuild) {
			if (debug) {
				printf("Qfunc_cgeneric: Rebuild Q-hash for id %d thread_id %d\n", id, thread_id);
			}
			if (a->Q[id]) {
				GMRFLib_free_tabulate_Qfunc(a->Q[id]);
			}
			double *a_tmp = Calloc(a->ntheta, double);
			for (int jj = 0; jj < a->ntheta; jj++) {
				a_tmp[jj] = a->theta[jj][thread_id][0];
				if (debug) {
					printf("\ttheta[%1d] %.12f\n", jj, a_tmp[jj]);
				}
			}

			if (debug) {
				printf("\tCall cgeneric\n");
			}
			x_out = a->model_func(INLA_CGENERIC_Q, a_tmp, a->data);
			if (a->debug) {
				inla_cgeneric_debug(stdout, a->secname, INLA_CGENERIC_Q, x_out);
			}

			assert(a->graph);
			assert((int) x_out[0] == -1);	       /* ONLY SUPPORT THIS, as its the same in any case... */
			if ((int) x_out[0] == -1) {
				// optimized output
				k = 1;
				len = (int) x_out[k++];
				assert(len == a->len_list);
				n = a->graph->n;
				GMRFLib_tabulate_Qfunc_from_list2(&(a->Q[id]), a->graph, a->len_list, a->ilist, a->jlist, &(x_out[k]), n, NULL);
			} else {
				k = 0;
				n = (int) x_out[k++];
				len = (int) x_out[k++];

				// we can overlay these arrays to avoid allocating new ones, since x_out is double
				ilist = (int *) &(x_out[k]);
				jlist = (int *) &(x_out[k + len]);
				Qijlist = (double *) &(x_out[k + 2 * len]);
				for (int jj = 0; jj < len; jj++) {
					ilist[jj] = (int) x_out[k + jj];
					jlist[jj] = (int) x_out[k + len + jj];
				}

				GMRFLib_tabulate_Qfunc_from_list2(&(a->Q[id]), a->graph, len, ilist, jlist, Qijlist, n, NULL);
				assert(a->graph->n == a->n);
			}
			Free(x_out);
			if (a->param[id]) {
				Memcpy(a->param[id], a_tmp, a->ntheta * sizeof(double));
				Free(a_tmp);
			} else {
				a->param[id] = a_tmp;
			}
			if (debug) {
				printf("\tRebuild for id %1d done\n", id);
			}
		}
	}

	if (j >= 0) {
		return (a->Q[id]->Qfunc(thread_id, i, j, NULL, a->Q[id]->Qfunc_arg));
	} else {
		return (a->Q[id]->Qfunc(thread_id, i, j, values, a->Q[id]->Qfunc_arg));
	}
}

double Qfunc_dmatern(int thread_id, int node, int nnode, double *UNUSED(values), void *arg)
{
	if (nnode < 0) {
		return NAN;
	}

	dmatern_arg_tp *a = (dmatern_arg_tp *) arg;
	double prec = map_exp(a->log_prec[thread_id][0], MAP_FORWARD, NULL);
	int rebuild, id = 0;
	const int debug = 0;

	GMRFLib_CACHE_SET_ID(id);
	rebuild = (a->param[id] == NULL || a->Q[id] == NULL);
	if (!rebuild) {
		// yes, log_prec is ...[0], so we start at 1
		rebuild = (a->param[id][1] != a->log_range[thread_id][0]) || (a->param[id][2] != a->log_nu[thread_id][0]);
	}

	if (rebuild) {
#pragma omp critical (Name_8b586dffc9258cff507f5f6098e240cbacd2c571)
		{
			// yes, log_prec is ...[0], so we start at 1
			double range, nu;
			if (debug) {
				printf("Qfunc_dmatern: Rebuild Q-hash for id %d\n", id);
			}

			if (!(a->Q[id])) {
				a->Q[id] = gsl_matrix_calloc(a->n, a->n);
			}

			a->param[id][1] = a->log_range[thread_id][0];
			a->param[id][2] = a->log_nu[thread_id][0];
			range = map_exp(a->param[id][1], MAP_FORWARD, NULL);
			nu = map_exp(a->param[id][2], MAP_FORWARD, NULL);

			if (debug) {
				printf("\trange %.4f nu %.4f\n", range, nu);
			}

			for (int i = 0; i < a->n; i++) {
				for (int j = i; j < a->n; j++) {
					double dist, val;

					dist = gsl_matrix_get(a->dist, i, j);
					val = inla_dmatern_cf(dist, range, nu);
					gsl_matrix_set(a->Q[id], i, j, val);
					gsl_matrix_set(a->Q[id], j, i, val);
				}
			}
			GMRFLib_gsl_spd_inverse(a->Q[id]);
		}
	}

	return prec * gsl_matrix_get(a->Q[id], node, nnode);
}

double mfunc_ar1(int thread_id, int UNUSED(i), void *arg)
{
	inla_ar1_arg_tp *a = (inla_ar1_arg_tp *) arg;
	return (a->mean[thread_id][0]);
}

double mfunc_rgeneric(int thread_id, int i, void *arg)
{
	inla_rgeneric_tp *a = (inla_rgeneric_tp *) arg;
	int rebuild, ii, id = 0;
	const int debug = 0;

	// possible fast return ?
	if (a->mu_zero) {
		return 0.0;
	}

	GMRFLib_CACHE_SET_ID(id);
	rebuild = (a->mu_param[id] == NULL || a->mu[id] == NULL);
	for (ii = 0; ii < a->ntheta && !rebuild; ii++) {
		rebuild = (a->mu_param[id][ii] != a->theta[ii][thread_id][0]);
	}

	if (rebuild) {
#pragma omp critical (Name_a878b76a6db370a6183df8d897b3a03b15039501)
		{
			int n, n_out, jj;
			double *x_out = NULL;

			if (debug) {
				printf("Rebuild mu-hash for id %d\n", id);
			}
			if (!(a->mu_param[id])) {
				a->mu_param[id] = Calloc(a->ntheta, double);
			}
			for (jj = 0; jj < a->ntheta; jj++) {
				a->mu_param[id][jj] = a->theta[jj][thread_id][0];
				if (debug) {
					printf("\ttheta[%1d] %.20g\n", jj, a->mu_param[id][jj]);
				}
			}

			if (debug) {
				printf("Call rgeneric\n");
			}
			inla_R_rgeneric(&n_out, &x_out, R_GENERIC_MU, a->model, &(a->ntheta), a->mu_param[id]);
			if (debug) {
				printf("Return from rgeneric with n_out= %1d\n", n_out);
			}

			assert(n_out > 0);
			n = (int) x_out[0];
			if (n > 0) {
				assert(n == a->n);
				if (!(a->mu[id])) {
					a->mu[id] = Calloc(n, double);
				}
				Memcpy(a->mu[id], &(x_out[1]), n * sizeof(double));
				a->mu_zero = 0;
			} else {
				a->mu_zero = 1;
			}
			Free(x_out);
		}
		// do a fast return here, so we do not need to allocate the a->mu[id] above. 
		if (a->mu_zero) {
			return 0.0;
		}
	}

	return (a->mu[id][i]);
}

double mfunc_cgeneric(int thread_id, int i, void *arg)
{
	inla_cgeneric_tp *a = (inla_cgeneric_tp *) arg;
	int rebuild, id = 0;
	const int debug = 0;

	// possible fast return ?
	if (a->mu_zero) {
		return 0.0;
	}

	GMRFLib_CACHE_SET_ID(id);
	rebuild = (a->mu_param[id] == NULL || a->mu[id] == NULL);
	for (int ii = 0; ii < a->ntheta && !rebuild; ii++) {
		rebuild = (a->mu_param[id][ii] != a->theta[ii][thread_id][0]);
	}

	if (rebuild) {
		int n;
		double *x_out = NULL;

		if (debug) {
			printf("Rebuild mu-hash for id %d thread_id %d\n", id, thread_id);
		}
		if (!(a->mu_param[id])) {
			assert(a->ntheta >= 0 && a->ntheta < 1000000);
			a->mu_param[id] = Calloc(a->ntheta, double);
		}
		for (int jj = 0; jj < a->ntheta; jj++) {
			a->mu_param[id][jj] = a->theta[jj][thread_id][0];
			if (debug) {
				printf("\ttheta[%1d] %.20g\n", jj, a->mu_param[id][jj]);
			}
		}

		if (debug) {
			printf("Call cgeneric in mfunc_cgeneric\n");
		}
		x_out = a->model_func(INLA_CGENERIC_MU, a->mu_param[id], a->data);
		if (debug) {
			printf("Return from cgeneric with x_out[0]= %1d\n", (int) x_out[0]);
		}

		n = (int) x_out[0];
		if (n > 0) {
			assert(n == a->n);
			if (!(a->mu[id])) {
				a->mu[id] = Calloc(n, double);
			}
			Memcpy(a->mu[id], &(x_out[1]), n * sizeof(double));
			a->mu_zero = 0;
		} else {
			a->mu_zero = 1;
		}
		Free(x_out);

		// do a fast return here, so we do not need to allocate the a->mu[id] above. 
		if (a->mu_zero) {
			return 0.0;
		}
	}

	return (a->mu[id][i]);
}

double Qfunc_clinear(int UNUSED(thread_id), int i, int j, double *UNUSED(values), void *arg)
{
	if (j < 0) {
		return NAN;
	}

	inla_clinear_tp *a = (inla_clinear_tp *) arg;
	assert(i == j);
	return (a->precision);
}

double Qfunc_sigm(int UNUSED(thread_id), int i, int j, double *UNUSED(values), void *arg)
{
	if (j < 0) {
		return NAN;
	}

	inla_sigm_tp *a = (inla_sigm_tp *) arg;
	assert(i == j);
	return (a->precision);
}

double Qfunc_log1exp(int UNUSED(thread_id), int i, int j, double *UNUSED(values), void *arg)
{
	if (j < 0) {
		return NAN;
	}

	inla_log1exp_tp *a = (inla_log1exp_tp *) arg;
	assert(i == j);
	return (a->precision);
}

double Qfunc_logdist(int UNUSED(thread_id), int i, int j, double *UNUSED(values), void *arg)
{
	if (j < 0) {
		return NAN;
	}

	inla_log1exp_tp *a = (inla_log1exp_tp *) arg;
	assert(i == j);
	return (a->precision);
}

double mfunc_clinear(int thread_id, int i, void *arg)
{
	inla_clinear_tp *a = (inla_clinear_tp *) arg;
	double beta = map_beta(a->beta[thread_id][0], MAP_FORWARD, a->beta_arg);

	return beta * a->x[i];
}

double mfunc_sigm(int thread_id, int i, void *arg)
{
	inla_sigm_tp *a = (inla_sigm_tp *) arg;
	double beta, halflife, shape, x, xx;

	beta = a->beta[thread_id][0];
	halflife = map_precision(a->log_halflife[thread_id][0], MAP_FORWARD, NULL);
	shape = map_precision(a->log_shape[thread_id][0], MAP_FORWARD, NULL);
	x = a->x[i];
	assert(x >= 0.0);

	xx = pow(x / halflife, shape);
	return beta * xx / (1.0 + xx);
}

double mfunc_revsigm(int thread_id, int i, void *arg)
{
	inla_sigm_tp *a = (inla_sigm_tp *) arg;
	double beta, halflife, shape, x, xx;

	beta = a->beta[thread_id][0];
	halflife = map_precision(a->log_halflife[thread_id][0], MAP_FORWARD, NULL);
	shape = map_precision(a->log_shape[thread_id][0], MAP_FORWARD, NULL);
	x = a->x[i];
	assert(x >= 0.0);

	xx = pow(x / halflife, shape);
	return beta * 1.0 / (1.0 + xx);
}

double mfunc_log1exp(int thread_id, int i, void *arg)
{
	inla_log1exp_tp *a = (inla_log1exp_tp *) arg;
	double beta, alpha, gama, x;

	beta = a->beta[thread_id][0];
	alpha = a->alpha[thread_id][0];
	gama = a->gamma[thread_id][0];
	x = a->x[i];
	assert(x >= 0.0);

	return beta * log1p(exp(alpha - gama * x));
}

double mfunc_logdist(int thread_id, int i, void *arg)
{
	inla_logdist_tp *a = (inla_logdist_tp *) arg;
	double beta, alpha1, alpha2, x;

	beta = a->beta[thread_id][0];
	alpha1 = map_exp(a->alpha1[thread_id][0], MAP_FORWARD, NULL);
	alpha2 = map_exp(a->alpha2[thread_id][0], MAP_FORWARD, NULL);
	x = a->x[i];
	assert(x > 0.0);

	return beta * (1.0 + exp(alpha1 * log(x) - alpha2 * x));
}

double Qfunc_mec(int thread_id, int i, int j, double *UNUSED(values), void *arg)
{
	if (j < 0) {
		return NAN;
	}

	inla_mec_tp *a = (inla_mec_tp *) arg;
	double prec_x = map_precision(a->log_prec_x[thread_id][0], MAP_FORWARD, NULL);
	double prec_obs = map_precision(a->log_prec_obs[thread_id][0], MAP_FORWARD, NULL);
	double beta = map_beta(a->beta[thread_id][0], MAP_FORWARD, a->map_beta_arg);
	double *scale = a->scale;

	assert(i == j);
	return (prec_x + prec_obs * scale[i]) / SQR(beta);
}

double mfunc_mec(int thread_id, int i, void *arg)
{
	inla_mec_tp *a = (inla_mec_tp *) arg;
	double prec_x = map_precision(a->log_prec_x[thread_id][0], MAP_FORWARD, NULL);
	double prec_obs = map_precision(a->log_prec_obs[thread_id][0], MAP_FORWARD, NULL);
	double beta = map_beta(a->beta[thread_id][0], MAP_FORWARD, a->map_beta_arg);
	double mean_x = map_identity(a->mean_x[thread_id][0], MAP_FORWARD, NULL);
	double *x_obs = a->x_obs;
	double *scale = a->scale;

	return beta * (prec_obs * scale[i] * x_obs[i] + prec_x * mean_x) / (prec_obs * scale[i] + prec_x);
}

double Qfunc_meb(int thread_id, int i, int j, double *UNUSED(values), void *arg)
{
	if (j < 0) {
		return NAN;
	}

	inla_meb_tp *a = (inla_meb_tp *) arg;

	double prec = map_precision(a->log_prec[thread_id][0], MAP_FORWARD, NULL);
	double beta = map_beta(a->beta[thread_id][0], MAP_FORWARD, a->map_beta_arg);
	double *scale = a->scale;

	assert(i == j);
	return prec * scale[i] / SQR(beta);
}

double mfunc_meb(int thread_id, int i, void *arg)
{
	inla_meb_tp *a = (inla_meb_tp *) arg;
	double beta = map_beta(a->beta[thread_id][0], MAP_FORWARD, a->map_beta_arg);

	return beta * a->x[i];
}

int inla_iid_wishart_nparam(int dim)
{
	/*
	 * return the number of theta parameters 
	 */
	return ((dim * (dim + 1)) / 2);
}

double Qfunc_iid_wishart(int thread_id, int node, int nnode, double *UNUSED(values), void *arg)
{
	if (nnode < 0) {
		return NAN;
	}

	/*
	 * This function returns the ij'th element of the precision matrix of the dim-dimensional iid. The parameterisation is given in the covariance matrix, so we
	 * need to compute the precision matrix. We store Q to avoid to compute it all the time.
	 */

	inla_iid_wishart_arg_tp *a = (inla_iid_wishart_arg_tp *) arg;
	int i, j, k, n_theta, dim, id = 0;
	const int debug = 0;
	double *vec = NULL;
	inla_wishart_hold_tp *hold = NULL;

	dim = a->dim;
	n_theta = inla_iid_wishart_nparam(a->dim);

	if (dim == 1) {
		/*
		 *  Fast return for the IID1D model; no need to do complicate things in this case
		 */
		return map_precision(a->log_prec[0][thread_id][0], MAP_FORWARD, NULL);
	}

	/*
	 * using this prevent us for using '#pragma omp critical' below, so its much quicker 
	 */
	GMRFLib_CACHE_SET_ID(id);
	assert(a->hold);
	hold = a->hold[id];
	if (hold == NULL) {
		a->hold[id] = Calloc(1, inla_wishart_hold_tp);
		a->hold[id]->vec = Calloc(n_theta, double);
		a->hold[id]->vec[0] = GMRFLib_uniform();
		a->hold[id]->Q = gsl_matrix_calloc(a->dim, a->dim);
		hold = a->hold[id];
	}

	vec = Calloc(n_theta, double);
	k = 0;
	for (i = 0; i < dim; i++) {
		vec[k] = map_precision(a->log_prec[i][thread_id][0], MAP_FORWARD, NULL);
		k++;
	}
	for (i = 0; i < n_theta - dim; i++) {
		vec[k] = map_rho(a->rho_intern[i][thread_id][0], MAP_FORWARD, NULL);
		k++;
	}
	assert(k == n_theta);

	if (memcmp((void *) vec, (void *) hold->vec, n_theta * sizeof(double))) {
		inla_iid_wishart_adjust(dim, vec);
		k = 0;
		for (i = 0; i < dim; i++) {
			gsl_matrix_set(hold->Q, i, i, 1.0 / vec[k]);
			k++;
		}
		for (i = 0; i < dim; i++) {
			for (j = i + 1; j < dim; j++) {
				double value = vec[k] / sqrt(vec[i] * vec[j]);
				gsl_matrix_set(hold->Q, i, j, value);
				gsl_matrix_set(hold->Q, j, i, value);
				k++;
			}
		}
		assert(k == n_theta);

		if (debug) {
			for (i = 0; i < n_theta; i++)
				printf("vec[%1d] = %.12f\n", i, vec[i]);
			GMRFLib_printf_gsl_matrix(stdout, hold->Q, " %.12f");
		}

		GMRFLib_gsl_spd_inverse(hold->Q);
		Memcpy((void *) hold->vec, (void *) vec, n_theta * sizeof(double));	/* YES! */
	}

	Free(vec);

	return gsl_matrix_get(hold->Q, node / a->n, nnode / a->n);
}

int inla_wishartk_build_Q(int dim, double *theta, gsl_matrix *Q, gsl_matrix *L)
{
	int i, j, k = 0, n_theta = INLA_WISHARTK_NTHETA(dim);
	const int debug = 0;
	gsl_matrix_set_zero(L);
	for (i = 0; i < dim; i++) {
		gsl_matrix_set(L, i, i, exp(theta[k++]));
	}
	for (j = 0; j < dim; j++) {
		for (i = j + 1; i < dim; i++) {
			gsl_matrix_set(L, i, j, theta[k++]);
		}
	}
	assert(k == n_theta);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, L, L, 0.0, Q);
	if (debug) {
		for (i = 0; i < n_theta; i++) {
			printf("theta[%1d] = %.8f\n", i, theta[i]);
		}
		FIXME("L");
		GMRFLib_printf_gsl_matrix(stdout, L, " %.6f");
		FIXME("Q");
		GMRFLib_printf_gsl_matrix(stdout, Q, " %.6f");
	}
	return GMRFLib_SUCCESS;
}

double Qfunc_iid_wishartk(int thread_id, int node, int nnode, double *UNUSED(values), void *arg)
{
	if (nnode < 0) {
		return NAN;
	}

	inla_iid_wishartk_arg_tp *a = (inla_iid_wishartk_arg_tp *) arg;
	int i, n_theta, dim, id = 0;
	inla_wishartk_hold_tp *hold = NULL;
	double *vec = NULL;

	dim = a->dim;
	n_theta = a->ntheta;

	GMRFLib_CACHE_SET_ID(id);
	hold = a->hold[id];
	if (hold == NULL) {
		a->hold[id] = Calloc(1, inla_wishartk_hold_tp);
		a->hold[id]->vec = Calloc(n_theta, double);
		a->hold[id]->vec[0] = GMRFLib_uniform();
		a->hold[id]->L = gsl_matrix_calloc(a->dim, a->dim);
		a->hold[id]->Q = gsl_matrix_calloc(a->dim, a->dim);
		hold = a->hold[id];
	}

	vec = a->vec[id];
	for (i = 0; i < n_theta; i++) {
		vec[i] = a->theta[i][thread_id][0];
	}

	if (memcmp((void *) vec, (void *) hold->vec, n_theta * sizeof(double))) {
		inla_wishartk_build_Q(dim, vec, hold->Q, hold->L);
		Memcpy((void *) hold->vec, (void *) vec, n_theta * sizeof(double));	/* YES! */
	}

	return gsl_matrix_get(hold->Q, node / a->n, nnode / a->n);
}

double Qfunc_iid2d(int thread_id, int i, int j, double *UNUSED(values), void *arg)
{
	if (j < 0) {
		return NAN;
	}

	inla_iid2d_arg_tp *a = (inla_iid2d_arg_tp *) arg;
	double prec0, prec1, rho;

	prec0 = map_precision(a->log_prec0[thread_id][0], MAP_FORWARD, NULL);
	prec1 = map_precision(a->log_prec1[thread_id][0], MAP_FORWARD, NULL);
	rho = map_rho(a->rho_intern[thread_id][0], MAP_FORWARD, NULL);

	if (i == j) {
		if (i < a->n) {
			return prec0 / (1.0 - SQR(rho));
		} else {
			return prec1 / (1.0 - SQR(rho));
		}
	}

	return -rho * sqrt(prec0 * prec1) / (1.0 - SQR(rho));
}

double Qfunc_2diid(int thread_id, int i, int j, double *UNUSED(values), void *arg)
{
	if (j < 0) {
		return NAN;
	}

	inla_2diid_arg_tp *a = (inla_2diid_arg_tp *) arg;
	double prec0, prec1, rho;

	prec0 = map_precision(a->log_prec0[thread_id][0], MAP_FORWARD, NULL);
	prec1 = map_precision(a->log_prec1[thread_id][0], MAP_FORWARD, NULL);
	rho = map_rho(a->rho_intern[thread_id][0], MAP_FORWARD, NULL);

	if (i == j) {
		if (GSL_IS_EVEN(i)) {
			return prec0 / (1.0 - SQR(rho));
		} else {
			return prec1 / (1.0 - SQR(rho));
		}
	}

	return -rho * sqrt(prec0 * prec1) / (1.0 - SQR(rho));
}

double Qfunc_2diid_wishart(int thread_id, int i, int j, double *UNUSED(values), void *arg)
{
	if (j < 0) {
		return NAN;
	}

	inla_2diid_arg_tp *a = (inla_2diid_arg_tp *) arg;
	double prec0, prec1, rho;

	prec0 = map_precision(a->log_prec0[thread_id][0], MAP_FORWARD, NULL);
	prec1 = map_precision(a->log_prec1[thread_id][0], MAP_FORWARD, NULL);
	rho = map_rho(a->rho_intern[thread_id][0], MAP_FORWARD, NULL);

	if (i == j) {
		if (i < a->n) {
			return prec0 / (1.0 - SQR(rho));
		} else {
			return prec1 / (1.0 - SQR(rho));
		}
	}

	return -rho * sqrt(prec0 * prec1) / (1.0 - SQR(rho));
}

double Qfunc_ar1(int thread_id, int i, int j, double *UNUSED(values), void *arg)
{
	if (j < 0) {
		return NAN;
	}

	inla_ar1_arg_tp *a = (inla_ar1_arg_tp *) arg;
	double phi, prec_marginal, prec, val;

	/*
	 * the log_prec is the log precision for the *marginal*; so we need to compute the log_prec for the innovation or conditional noise. 
	 */
	phi = map_phi(a->phi_intern[thread_id][0], MAP_FORWARD, NULL);
	prec_marginal = map_precision(a->log_prec[thread_id][0], MAP_FORWARD, NULL);
	prec = prec_marginal / (1.0 - SQR(phi));

	if (a->cyclic) {
		if (i == j) {
			val = prec * (1.0 + SQR(phi));
		} else {
			val = -prec * phi;
		}
	} else {
		if (i != j) {
			val = -prec * phi;
		} else {
			if (i == 0 || i == a->n - 1) {
				val = prec;
			} else {
				val = prec * (1.0 + SQR(phi));
			}
		}
	}

	return val;
}

double Qfunc_ar1c(int thread_id, int i, int j, double *UNUSED(values), void *arg)
{
	if (j < 0) {
		return NAN;
	}

	inla_ar1c_arg_tp *a = (inla_ar1c_arg_tp *) arg;
	double phi, prec_marginal, prec, val;
	int ii, jj, m;

	/*
	 * the log_prec is the log precision for the *marginal*; so we need to compute the log_prec for the innovations
	 */
	phi = map_phi(a->phi_intern[thread_id][0], MAP_FORWARD, NULL);
	prec_marginal = map_precision(a->log_prec[thread_id][0], MAP_FORWARD, NULL);
	prec = prec_marginal / (1.0 - SQR(phi));

	ii = IMIN(i, j);
	jj = IMAX(i, j);

	if (jj < a->n) {
		// the AR1 part
		if (ii != jj) {
			val = -prec * phi;
		} else {
			if (ii == 0 || ii == a->n - 1) {
				val = prec;
			} else {
				val = prec * (1.0 + SQR(phi));
			}
		}
	} else if (ii < a->n) {
		// the interaction
		m = jj - a->n;
		if (ii == 0) {
			val = prec * phi * GMRFLib_matrix_get(ii, m, a->Z);
		} else if (ii == a->n - 1) {
			val = -prec * GMRFLib_matrix_get(ii - 1, m, a->Z);
		} else {
			val = prec * (phi * GMRFLib_matrix_get(ii, m, a->Z) - GMRFLib_matrix_get(ii - 1, m, a->Z));
		}
	} else {
		// the beta-block
		int iii = ii - a->n, jjj = jj - a->n;
		val = GMRFLib_matrix_get(iii, jjj, a->Qbeta) + prec * GMRFLib_matrix_get(iii, jjj, a->ZZ);
	}

	return val;
}

double Qfunc_ou(int thread_id, int i, int j, double *UNUSED(values), void *arg)
{
	if (j < 0) {
		return NAN;
	}

	inla_ou_arg_tp *a = (inla_ou_arg_tp *) arg;
	double phi, prec, w, v, delta;

	phi = map_exp(a->phi_intern[thread_id][0], MAP_FORWARD, NULL);
	prec = map_precision(a->log_prec[thread_id][0], MAP_FORWARD, NULL);

	if (i != j) {
		int ii = IMAX(i, j);
		delta = a->locations[ii] - a->locations[ii - 1];
		w = 1.0 / ONE_mexp(-2.0 * phi * delta);
		v = exp(-phi * delta);

		return -prec * w * v;
	} else {
		if (i == 0) {
			delta = a->locations[i + 1] - a->locations[i];
			w = 1.0 / ONE_mexp(-2.0 * phi * delta);
			v = exp(-phi * delta);

			return prec * (1.0 + w * SQR(v));
		} else if (i == a->n - 1) {
			delta = a->locations[i] - a->locations[i - 1];
			w = 1.0 / ONE_mexp(-2.0 * phi * delta);

			return prec * w;
		} else {
			delta = a->locations[i + 1] - a->locations[i];
			w = 1.0 / ONE_mexp(-2.0 * phi * delta);
			v = exp(-phi * delta);

			double ddelta = a->locations[i] - a->locations[i - 1];
			double ww = 1.0 / ONE_mexp(-2.0 * phi * ddelta);

			return prec * (ww + w * SQR(v));
		}
	}
	abort();
	return 0.0;
}

double Qfunc_besag(int thread_id, int i, int j, double *UNUSED(values), void *arg)
{
	if (j < 0) {
		return NAN;
	}

	inla_besag_Qfunc_arg_tp *a;
	double prec;

	a = (inla_besag_Qfunc_arg_tp *) arg;
	prec = (a->log_prec ? map_precision(a->log_prec[thread_id][0], MAP_FORWARD, NULL) : 1.0);

	if (a->prec_scale) {
		if (a->prec_scale[i] > 0.0) {
			prec *= a->prec_scale[i];
			// normal return
			return prec * (i == j ? a->graph->nnbs[i] : -1.0);
		} else if (a->prec_scale[i] < 0.0) {
			// node with no neibours and marg variance = 1
			assert((i == j) && a->graph->nnbs[i] == 0);
			return (prec);
		} else if (a->prec_scale[i] == 0.0) {
			// uniform
			return 0.0;
		} else {
			assert(0 == 1);
		}
	}
	// ``classical model''
	return prec * (i == j ? a->graph->nnbs[i] : -1.0);
}

double Qfunc_besag2(int thread_id, int i, int j, double *UNUSED(values), void *arg)
{
	if (j < 0) {
		return NAN;
	}

	inla_besag2_Qfunc_arg_tp *aa;
	double a;

	aa = (inla_besag2_Qfunc_arg_tp *) arg;
	if (aa->log_a) {
		a = map_exp(aa->log_a[thread_id][0], MAP_FORWARD, NULL);
	} else {
		a = 1.0;
	}

	if (IMAX(i, j) < aa->besag_arg->graph->n) {
		if (i == j) {
			return (Qfunc_besag(thread_id, i, j, NULL, (void *) aa->besag_arg) + aa->precision / SQR(a)) / SQR(a) + aa->diag;
		} else {
			return Qfunc_besag(thread_id, i, j, NULL, (void *) aa->besag_arg) / SQR(a);
		}
	} else if (IMIN(i, j) >= aa->besag_arg->graph->n) {
		return aa->precision;
	} else {
		return -aa->precision / SQR(a);
	}
}

double Qfunc_besagproper(int thread_id, int i, int j, double *UNUSED(values), void *arg)
{
	if (j < 0) {
		return NAN;
	}

	inla_besag_proper_Qfunc_arg_tp *a;
	double prec;

	a = (inla_besag_proper_Qfunc_arg_tp *) arg;
	prec = map_precision(a->log_prec[thread_id][0], MAP_FORWARD, NULL);
	if (i == j) {
		double diag = map_exp(a->log_diag[thread_id][0], MAP_FORWARD, NULL);
		return prec * (diag + a->graph->nnbs[i]);
	} else {
		return -prec;
	}
}

double Qfunc_besagproper2(int thread_id, int i, int j, double *UNUSED(values), void *arg)
{
	if (j < 0) {
		return NAN;
	}

	inla_besag_proper2_Qfunc_arg_tp *a;
	double prec;
	double lambda;

	a = (inla_besag_proper2_Qfunc_arg_tp *) arg;
	prec = map_precision(a->log_prec[thread_id][0], MAP_FORWARD, NULL);
	lambda = map_probability(a->logit_lambda[thread_id][0], MAP_FORWARD, NULL);
	if (i == j) {
		return prec * ((1.0 - lambda) + lambda * a->graph->nnbs[i]);
	} else {
		return -prec * lambda;
	}
}

double Qfunc_intslope(int thread_id, int node, int nnode, double *UNUSED(values), void *arg)
{
	if (nnode < 0) {
		return NAN;
	}

	int i, imin, imax, idx, subject, strata, icase;
	double val = 0.0, xval = 0.0, gam, z = NAN;
	inla_intslope_arg_tp *a = (inla_intslope_arg_tp *) arg;

	imin = IMIN(node, nnode);
	imax = IMAX(node, nnode);

	if (imax < a->n) {
		// the diagonal part
		return a->precision;
	} else if (imin >= a->n) {
		// the Wishart part
		imin -= a->n;
		imax -= a->n;
		val = Qfunc_iid_wishart(thread_id, imin, imax, NULL, (void *) a->warg);
		subject = (imin == imax ? (imin < a->nsubject ? imin : imin - a->nsubject) : imin);
		icase = (imax < a->nsubject ? 0 : (imin >= a->nsubject ? 1 : 2));
		for (i = 0; i < a->subject_idx[subject]->n; i++) {
			idx = a->subject_idx[subject]->idx[i];
			strata = (int) GMRFLib_matrix_get(idx, INTSLOPE_STRATA, a->def);
			gam = a->theta_gamma[strata][thread_id][0];
			if (icase > 0) {
				z = GMRFLib_matrix_get(idx, INTSLOPE_Z, a->def);
			}
			switch (icase) {
			case 0:
			{
				xval += SQR(gam);
			}
				break;
			case 1:
			{
				xval += SQR(gam * z);
			}
				break;
			case 2:
			{
				xval += SQR(gam) * z;
			}
				break;
			}
		}

		return val + a->precision * xval;
	} else {
		imax -= a->n;
		strata = (int) GMRFLib_matrix_get(imin, INTSLOPE_STRATA, a->def);
		gam = a->theta_gamma[strata][thread_id][0];
		if (imax < a->nsubject) {
			val = -gam;
		} else {
			z = GMRFLib_matrix_get(imin, INTSLOPE_Z, a->def);
			val = -gam * z;
		}
		return a->precision * val;
	}

	assert(0 == 1);
	return 0.0;
}

double iid_mfunc(int idx, void *UNUSED(arg))
{
	return 1.0 + idx;
}

double Qfunc_copy_part00(int thread_id, int i, int j, double *UNUSED(values), void *arg)
{
	if (j < 0) {
		return NAN;
	}

	inla_copy_arg_tp *a = (inla_copy_arg_tp *) arg;

	if (i == j) {
		double beta = a->map_beta(a->beta[thread_id][0], MAP_FORWARD, a->map_beta_arg);
		return a->Qfunc(thread_id, i, j, NULL, a->Qfunc_arg) + a->precision * SQR(beta);
	} else {
		return a->Qfunc(thread_id, i, j, NULL, a->Qfunc_arg);
	}
}

double Qfunc_copy_part01(int thread_id, int UNUSED(i), int j, double *UNUSED(values), void *arg)
{
	if (j < 0) {
		return NAN;
	}

	inla_copy_arg_tp *a = (inla_copy_arg_tp *) arg;
	double beta = a->map_beta(a->beta[thread_id][0], MAP_FORWARD, a->map_beta_arg);

	return -a->precision * beta;
}

double Qfunc_copy_part11(int UNUSED(thread_id), int UNUSED(i), int j, double *UNUSED(values), void *arg)
{
	if (j < 0) {
		return NAN;
	}

	inla_copy_arg_tp *a = (inla_copy_arg_tp *) arg;

	return a->precision;
}

double Qfunc_scopy_part00(int thread_id, int i, int j, double *UNUSED(values), void *arg)
{
	if (j < 0) {
		return NAN;
	}

	inla_scopy_arg_tp *a = (inla_scopy_arg_tp *) arg;

	int cache_idx = 0;
	GMRFLib_CACHE_SET_ID(cache_idx);

	int build = 0;
	for (int k = 0; k < a->nbeta; k++) {
		if (a->betas[k][thread_id][0] != a->cache00[cache_idx]->betas[k]) {
			build = 1;
			break;
		}
	}

	if (build) {
#pragma omp critical (Name_a1fd55509a70889a1417fb08664e246f4517ee2a)
		{
			assert(a->nbeta == a->W->ncol);
			for (int k = 0; k < a->nbeta; k++) {
				a->cache00[cache_idx]->betas_tmp[k] = a->betas[k][thread_id][0];
			}

			if (a->nbeta > 2) {
				double *theta = Calloc(a->nbeta, double);
				for (int k = 0; k < a->nbeta; k++) {
					double *b = a->cache00[cache_idx]->betas_tmp;
					for (int jj = 0; jj < a->nbeta; jj++) {
						theta[k] += a->W->A[k + jj * a->nbeta] * b[jj];
					}
				}
				GMRFLib_spline_free(a->cache00[cache_idx]->splinefun);
				a->cache00[cache_idx]->splinefun = GMRFLib_spline_create(a->loc_beta, theta, a->nbeta);
				Free(theta);
			}
			Memcpy(a->cache00[cache_idx]->betas, a->cache00[cache_idx]->betas_tmp, a->nbeta * sizeof(double));
		}
	}

	if (i == j) {
		if (a->nbeta == 2) {
			double *ab = a->cache00[cache_idx]->betas;
			double beta = ab[0] + ab[1] * (a->cov_beta[i] - a->loc_mid) / a->loc_len;
			return a->Qfunc(thread_id, i, j, NULL, a->Qfunc_arg) + a->precision * SQR(beta);
		} else {
			double beta = GMRFLib_spline_eval(a->cov_beta[i], a->cache00[cache_idx]->splinefun);
			return a->Qfunc(thread_id, i, j, NULL, a->Qfunc_arg) + a->precision * SQR(beta);
		}
	} else {
		return a->Qfunc(thread_id, i, j, NULL, a->Qfunc_arg);
	}
}

double Qfunc_scopy_part01(int thread_id, int i, int j, double *UNUSED(values), void *arg)
{
	if (j < 0) {
		return NAN;
	}

	inla_scopy_arg_tp *a = (inla_scopy_arg_tp *) arg;

	int cache_idx = 0;
	GMRFLib_CACHE_SET_ID(cache_idx);

	int build = 0;
	for (int k = 0; k < a->nbeta; k++) {
		if (a->betas[k][thread_id][0] != a->cache01[cache_idx]->betas[k]) {
			build = 1;
			break;
		}
	}

	if (build) {
#pragma omp critical (Name_562559af23fb070f255b55089f79f0c69a8b73a2)
		{
			assert(a->nbeta == a->W->ncol);
			for (int k = 0; k < a->nbeta; k++) {
				a->cache01[cache_idx]->betas_tmp[k] = a->betas[k][thread_id][0];
			}

			if (a->nbeta > 2) {
				double *theta = Calloc(a->nbeta, double);
				for (int k = 0; k < a->nbeta; k++) {
					double *b = a->cache01[cache_idx]->betas_tmp;
					for (int jj = 0; jj < a->nbeta; jj++) {
						theta[k] += a->W->A[k + jj * a->nbeta] * b[jj];
					}
				}
				GMRFLib_spline_free(a->cache01[cache_idx]->splinefun);
				a->cache01[cache_idx]->splinefun = GMRFLib_spline_create(a->loc_beta, theta, a->nbeta);
				Free(theta);
			}
			Memcpy(a->cache01[cache_idx]->betas, a->cache01[cache_idx]->betas_tmp, a->nbeta * sizeof(double));
		}
	}

	assert(i == j);

	double beta_i;
	if (a->nbeta == 2) {
		double *ab = a->cache01[cache_idx]->betas;
		beta_i = ab[0] + ab[1] * (a->cov_beta[i] - a->loc_mid) / a->loc_len;
	} else {
		beta_i = GMRFLib_spline_eval(a->cov_beta[i], a->cache01[cache_idx]->splinefun);
	}
	return -a->precision * beta_i;
}

double Qfunc_scopy_part11(int UNUSED(thread_id), int UNUSED(i), int j, double *UNUSED(values), void *arg)
{
	if (j < 0) {
		return NAN;
	}
	inla_scopy_arg_tp *a = (inla_scopy_arg_tp *) arg;
	return a->precision;
}

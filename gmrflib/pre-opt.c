
/* pre-opt.c
 * 
 * Copyright (C) 2021-2021 Havard Rue
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
#ifndef GITCOMMIT
#define GITCOMMIT
#endif
static const char GitID[] = "file: " __FILE__ "  " GITCOMMIT;

#include <time.h>
#include <strings.h>
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif
#include <stdlib.h>
#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"


int GMRFLib_preopt_init(GMRFLib_preopt_tp **preopt, 
			int nf, 
			GMRFLib_graph_tp ** f_graph, GMRFLib_Qfunc_tp ** f_Qfunc,
			void **f_Qfunc_arg, char *f_sumzero, GMRFLib_constr_tp ** f_constr,
			GMRFLib_Qfunc_tp *** ff_Qfunc, void ***ff_Qfunc_arg,
			int nbeta, double **covariate, double *prior_precision, 
			GMRFLib_ai_param_tp * UNUSED(ai_par))
{
	int i, ii, j, k, N, *idx_map_f = NULL, *idx_map_beta = NULL, offset; 
	GMRFLib_preopt_arg_tp *arg = NULL;
	GMRFLib_constr_tp *fc = NULL;

	if (!preopt) {
		return GMRFLib_SUCCESS;
	}

	*preopt = Calloc(1, GMRFLib_preopt_tp);
	arg = Calloc(1, GMRFLib_preopt_arg_tp);
	nbeta = IMAX(0, nbeta);
	nf = IMAX(0, nf);
	arg->nf = nf;
	arg->f_Qfunc = f_Qfunc;
	arg->f_Qfunc_arg = f_Qfunc_arg;
	arg->ff_Qfunc = ff_Qfunc;
	arg->ff_Qfunc_arg = ff_Qfunc_arg;
	arg->f_graph = f_graph;
	arg->nbeta = nbeta;
	arg->covariate = covariate;
	arg->prior_precision = prior_precision;

	if (ff_Qfunc) {
		/*
		 * check that the specification is symmetric, as the implementation depends on it. 
		 */
		for (i = 0; i < nf; i++) {
			for (j = 0; j < nf; j++) {
				if (i != j) {
					if (ff_Qfunc[i][j]) {
						GMRFLib_ASSERT(ff_Qfunc[i][j] == ff_Qfunc[j][i], GMRFLib_EPARAMETER);
					}
					if (ff_Qfunc_arg) {
						GMRFLib_ASSERT(ff_Qfunc_arg[i][j] == ff_Qfunc_arg[j][i], GMRFLib_EPARAMETER);
					}
				}
			}
		}
	}


	/*
	 * Our first job, is to go through the model and compute all interactions etc that are defined through the \eta-model. 
	 * define the index-mapping. The outline of x is (eta, f[0], ..., beta[0], ...)
	 */
	offset = 0;
	if (nf) {
		idx_map_f = Calloc(nf + 1, int);
		for (i = 0; i < nf; i++) {
			idx_map_f[i] = offset;
			offset += f_graph[i]->n;
		}
		idx_map_f[nf] = offset;
	}
	if (nbeta) {
		idx_map_beta = Calloc(nbeta + 1, int);
		for (i = 0; i < nbeta; i++) {
			idx_map_beta[i] = offset;
			offset++;
		}
		idx_map_beta[nbeta] = offset;
	}
	N = offset;					       /* N is the grand-total. */
	P(N);
	
	/*
	 * If we have cross-terms, make sure to mark these cross-terms as neigbours; just fill them with zero's.
	 * If they are connected in the data, then this value
	 * will be overrided and this is how it should be.
	 */
	GMRFLib_ged_tp *ged = NULL;
	GMRFLib_ged_init(&ged, NULL);
	GMRFLib_ged_add(ged, N-1, N-1);

	if (nf) {
		for (j = 0; j < nf; j++) {
			GMRFLib_ged_insert_graph(ged, f_graph[j], idx_map_f[j]);
		}
	}
	for (i = 0; i < nf; i++) {
		for (j = i + 1; j < nf; j++) {
			if (ff_Qfunc && ff_Qfunc[i][j]) {
				for (k = 0; k < IMIN(f_graph[i]->n, f_graph[j]->n); k++) {
					GMRFLib_ged_add(ged, idx_map_f[i] + k, idx_map_f[j] + k);
				}
			}
		}
	}

	if (nbeta) {
		for(i = 0; i < nbeta; i++) {
			GMRFLib_ged_add(ged, idx_map_beta[i], idx_map_beta[i]);
		}
	}

	GMRFLib_ged_build(&((*preopt)->latent_graph), ged);
	GMRFLib_ged_free(ged);

	(*preopt)->latent_Qfunc = GMRFLib_preopt_Qfunc;
	(*preopt)->latent_Qfunc_arg = (void *) arg;

	/*
	 * build the constraint, if any. Only simple sum-to-zero constraints are supported.
	 */
	int nconstr = 0;
	if (nf && f_sumzero) {
		for (k = 0; k < nf; k++) {
			nconstr += (f_sumzero[k] ? 1 : 0);
		}
	}

	if (nf && f_constr) {
		for (k = 0; k < nf; k++) {
			fc = f_constr[k];
			if (fc) {
				nconstr += fc->nc;
			}
		}
	}
	if (nconstr) {
		GMRFLib_constr_tp *constr = NULL;
		int constr_no;
		GMRFLib_make_empty_constr(&constr);
		constr->a_matrix = Calloc((*preopt)->latent_graph->n * nconstr, double);
		constr->e_vector = Calloc(nconstr, double);
		constr->nc = nconstr;
		constr_no = 0;
		if (nf && f_sumzero) {
			for (k = 0; k < nf; k++) {
				if (f_sumzero[k]) {
					for (i = idx_map_f[k]; i < idx_map_f[k + 1]; i++) {
						constr->a_matrix[i * constr->nc + constr_no] = 1.0;
					}
					constr->e_vector[constr_no] = 0.0;
					constr_no++;
				}
			}
		}
		if (nf && f_constr) {
			for (k = 0; k < nf; k++) {
				fc = f_constr[k];
				if (fc) {
					for (j = 0; j < fc->nc; j++) {
						for (i = idx_map_f[k], ii = 0; i < idx_map_f[k + 1]; i++, ii++) {
							constr->a_matrix[i * constr->nc + constr_no] = fc->a_matrix[ii * fc->nc + j];
						}
						constr->e_vector[constr_no] = fc->e_vector[j];
						constr_no++;
					}
				}
			}
		}

		GMRFLib_ASSERT(constr_no == constr->nc, GMRFLib_ESNH);
		GMRFLib_prepare_constr(constr, (*preopt)->latent_graph, 0);
		(*preopt)->latent_constr = constr;
	} else {
		(*preopt)->latent_constr = NULL;
	}

	arg->idx_map_f = idx_map_f;
	arg->idx_map_beta = idx_map_beta;
	arg->n = (*preopt)->latent_graph->n;
	arg->what_type = Calloc(N, GMRFLib_preopt_type_tp);
	for (i = 0; i < N; i++) {
		arg->what_type[i] = GMRFLib_preopt_what_type(i, arg);
	}

	GMRFLib_openmp_implement_strategy(GMRFLib_OPENMP_PLACES_DEFAULT, NULL, NULL);

	return GMRFLib_SUCCESS;
}

GMRFLib_preopt_type_tp GMRFLib_preopt_what_type(int node, GMRFLib_preopt_arg_tp * a)
{
	int i;
	GMRFLib_preopt_type_tp t = { GMRFLib_PREOPT_TP___VOID, -1, -1 };

	if (a->nf && node < a->idx_map_f[a->nf]) {
		t.tp = GMRFLib_PREOPT_TP_F;
		for (i = 0; i < a->nf; i++) {
			if (node < a->idx_map_f[i + 1]) {
				t.tp_idx = i;
				t.idx = node - a->idx_map_f[i];
				break;
			}
		}
	} else if (a->nbeta && node < a->idx_map_beta[a->nbeta]) {
		t.tp = GMRFLib_PREOPT_TP_BETA;
		t.tp_idx = node - a->idx_map_beta[0];
		t.idx = 0;
	} else {
		GMRFLib_ASSERT_RETVAL(0 == 1, GMRFLib_ESNH, t);
	}

	return t;
}

double GMRFLib_preopt_Qfunc(int node, int nnode, double *UNUSED(values), void *arg)
{
	if (node >= 0 && nnode < 0) {
		return NAN;
	}

	/*
	 * this is Qfunction for the preopt-function 
	 */
	GMRFLib_preopt_arg_tp *a = NULL;
	int ii, jj;
	double value = 0.0;
	GMRFLib_preopt_type_tp it, jt;

	a = (GMRFLib_preopt_arg_tp *) arg;
	ii = IMIN(node, nnode);
	jj = IMAX(node, nnode);

	// old non-caching code: it = GMRFLib_preopt_what_type(ii, a); jt = GMRFLib_preopt_what_type(jj, a);
	it = a->what_type[ii];
	jt = a->what_type[jj];

	switch (it.tp) {
	case GMRFLib_PREOPT_TP_F:
		switch (jt.tp) {
		case GMRFLib_PREOPT_TP_F:
			if (it.tp_idx == jt.tp_idx) {
				if ((it.idx == jt.idx) || GMRFLib_graph_is_nb(it.idx, jt.idx, a->f_graph[it.tp_idx])) {
					value += a->f_Qfunc[it.tp_idx] (it.idx, jt.idx, NULL, (a->f_Qfunc_arg ? a->f_Qfunc_arg[it.tp_idx] : NULL));
				}
			}
			/*
			 * only for the same index and different types; used to define `interaction between fields'. this is a 'workaround' for a INLA problem.. 
			 */
			if (a->ff_Qfunc) {
				if ((it.idx == jt.idx) && (it.tp_idx != jt.tp_idx) && a->ff_Qfunc[it.tp_idx][jt.tp_idx]) {
					value +=
						a->ff_Qfunc[it.tp_idx][jt.tp_idx] (it.idx, jt.idx, NULL,
										   (a->ff_Qfunc_arg ? a->ff_Qfunc_arg[it.tp_idx][jt.tp_idx] : NULL));
				}
			}
			return value;
		case GMRFLib_PREOPT_TP_BETA:
			return value;
		default:
			GMRFLib_ASSERT_RETVAL(0 == 1, GMRFLib_ESNH, 0.0);
		}
		GMRFLib_ASSERT_RETVAL(0 == 1, GMRFLib_ESNH, 0.0);
		break;

	case GMRFLib_PREOPT_TP_BETA:
		switch (jt.tp) {
		case GMRFLib_PREOPT_TP_BETA:
			if (it.tp_idx == jt.tp_idx) {
				value += (a->prior_precision ? a->prior_precision[it.tp_idx] : 0.0);
			}
			return value;
		default:
			GMRFLib_ASSERT_RETVAL(0 == 1, GMRFLib_ESNH, 0.0);
		}
		GMRFLib_ASSERT_RETVAL(0 == 1, GMRFLib_ESNH, 0.0);
		break;

	default:
		GMRFLib_ASSERT_RETVAL(0 == 1, GMRFLib_ESNH, 0.0);
		break;
	}

	return value;
}

int GMRFLib_free_preopt(GMRFLib_preopt_tp * preopt)
{
	if (!preopt) {
		return GMRFLib_SUCCESS;
	}
	GMRFLib_preopt_arg_tp *a = (GMRFLib_preopt_arg_tp *) preopt->latent_Qfunc_arg;

	GMRFLib_graph_free(preopt->latent_graph);
	GMRFLib_free_constr(preopt->latent_constr);

	Free(a->idx_map_f);
	Free(a->idx_map_beta);
	Free(a->what_type);

	Free(a);
	Free(preopt);

	return GMRFLib_SUCCESS;
}

int GMRFLib_preopt_test(GMRFLib_preopt_tp *preopt) 
{
	if (!preopt) {
		return GMRFLib_SUCCESS;
	}

	GMRFLib_printf_graph(stdout, preopt->latent_graph);
	GMRFLib_printf_Qfunc(stdout, preopt->latent_graph, preopt->latent_Qfunc, preopt->latent_Qfunc_arg);
	GMRFLib_printf_constr(stdout, preopt->latent_constr, preopt->latent_graph); 
	return GMRFLib_SUCCESS;
}


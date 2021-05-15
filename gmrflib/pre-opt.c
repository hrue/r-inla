
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


int GMRFLib_preopt_init(GMRFLib_preopt_tp ** preopt,
			int nlike, int nf, int **c, double **w,
			GMRFLib_graph_tp ** f_graph, GMRFLib_Qfunc_tp ** f_Qfunc,
			void **f_Qfunc_arg, char *f_sumzero, GMRFLib_constr_tp ** f_constr,
			GMRFLib_Qfunc_tp *** ff_Qfunc, void ***ff_Qfunc_arg,
			int nbeta, double **covariate, double *prior_precision, GMRFLib_bfunc_tp ** bfunc, GMRFLib_ai_param_tp * UNUSED(ai_par))
{
	int i, ii, j, jj, k, N, *idx_map_f = NULL, *idx_map_beta = NULL, offset;
	int debug = 0;
	GMRFLib_preopt_arg_tp *arg = NULL;
	GMRFLib_constr_tp *fc = NULL;

	if (!preopt) {
		return GMRFLib_SUCCESS;
	}

	double **ww = NULL;
	ww = Calloc(nf, double *);
	for (i = 0; i < nf; i++) {
		ww[i] = Calloc(nlike, double);
		for (j = 0; j < nlike; j++) {
			ww[i][j] = 1.0;
		}
	}
	if (w) {
		for (i = 0; i < nf; i++) {
			if (w[i]) {
				memcpy(ww[i], w[i], nlike * sizeof(double));
			}
		}
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

	/*
	 * If we have cross-terms, make sure to mark these cross-terms as neigbours; just fill them with zero's.
	 * If they are connected in the data, then this value
	 * will be overrided and this is how it should be.
	 */
	GMRFLib_ged_tp *ged = NULL;
	GMRFLib_ged_init(&ged, NULL);
	GMRFLib_ged_add(ged, N - 1, N - 1);

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
		for (i = 0; i < nbeta; i++) {
			GMRFLib_ged_add(ged, idx_map_beta[i], idx_map_beta[i]);
		}
	}

	GMRFLib_ged_build(&((*preopt)->latent_graph), ged);
	GMRFLib_ged_free(ged);

	arg->latent_graph = (*preopt)->latent_graph;	       /* just a copy */
	(*preopt)->latent_Qfunc = GMRFLib_preopt_latent_Qfunc;
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
	arg->what_type = Calloc(N, GMRFLib_preopt_type_tp);
	for (i = 0; i < N; i++) {
		arg->what_type[i] = GMRFLib_preopt_what_type(i, arg);
	}


	if (debug) {
		printf("\tndata %1d nf %1d nbeta %1d\n", nlike, nf, nbeta);
		for (i = 0; i < nlike; i++) {
			printf("data %1d\n", i);
			for (j = 0; j < nf; j++) {
				printf("\t\tf[%1d]  index %1d  weight %.6f\n", j, c[j][i], ww[j][i]);
			}
			for (j = 0; j < nbeta; j++) {
				printf("\t\tbeta[%1d]  x %.6f\n", j, covariate[j][i]);
			}
		}
	}
	// build up structure for the likelihood part

	GMRFLib_idxval_tp **At_idxval = NULL;
	At_idxval = Calloc(N, GMRFLib_idxval_tp *);
	arg->At_idxval = At_idxval;

	for (i = 0; i < nlike; i++) {
		int idx;
		double val;

		for (jj = 0; jj < nf; jj++) {
			if (c[jj][i] >= 0 && ww[jj][i]) {
				idx = c[jj][i] + idx_map_f[jj];
				val = ww[jj][i];
				GMRFLib_idxval_add(&(At_idxval[idx]), i, val);
			}
		}
		for (jj = 0; jj < nbeta; jj++) {
			if (covariate[jj][i]) {
				idx = idx_map_beta[jj];
				val = covariate[jj][i];
				GMRFLib_idxval_add(&(At_idxval[idx]), i, val);
			}
		}
	}

	ged = NULL;
	GMRFLib_ged_init(&ged, NULL);
	for (i = 0; i < N; i++) {
		GMRFLib_ged_add(ged, i, i);
	}

	for (i = 0; i < nlike; i++) {
		GMRFLib_idx_tp *idx = NULL;
		for (jj = 0; jj < nf; jj++) {
			if (c[jj][i] >= 0 && ww[jj][i]) {
				GMRFLib_idx_add(&idx, c[jj][i] + idx_map_f[jj]);
			}
		}
		for (jj = 0; jj < nbeta; jj++) {
			if (covariate[jj][i]) {
				GMRFLib_idx_add(&idx, idx_map_beta[jj]);
			}
		}
		for (j = 0; j < idx->n; j++) {
			for (jj = j + 1; jj < idx->n; jj++) {
				GMRFLib_ged_add(ged, idx->idx[j], idx->idx[jj]);
			}
		}
		GMRFLib_idx_free(idx);
	}

	GMRFLib_graph_tp *g = NULL;
	GMRFLib_ged_build(&g, ged);
	GMRFLib_ged_free(ged);
	assert(g->n == N);

	GMRFLib_idxval_tp ***AtA_idxval = NULL;
	AtA_idxval = Calloc(N, GMRFLib_idxval_tp **);

	for (i = 0; i < g->n; i++) {
		AtA_idxval[i] = Calloc(1 + g->lnnbs[i], GMRFLib_idxval_tp *);
		GMRFLib_idxval_add(&(AtA_idxval[i][0]), i, 0.0);
	}

	for (i = 0; i < nlike; i++) {
		GMRFLib_idxval_tp *idxval = NULL;

		for (jj = 0; jj < nf; jj++) {
			if (c[jj][i] >= 0 && ww[jj][i]) {
				GMRFLib_idxval_add(&idxval, c[jj][i] + idx_map_f[jj], ww[jj][i]);
			}
		}
		for (jj = 0; jj < nbeta; jj++) {
			if (covariate[jj][i]) {
				GMRFLib_idxval_add(&idxval, idx_map_beta[jj], covariate[jj][i]);
			}
		}

		for (j = 0; j < idxval->n; j++) {
			for (jj = j; jj < idxval->n; jj++) {
				int imin, imax, index;
				imin = IMIN(idxval->store[j].idx, idxval->store[jj].idx);
				imax = IMAX(idxval->store[j].idx, idxval->store[jj].idx);

				if (imin == imax) {
					index = 0;
				} else {
					index = 1 + GMRFLib_iwhich_sorted(imax, g->lnbs[imin], g->lnnbs[imin]);
					assert(index > 0);
				}
				GMRFLib_idxval_add(&(AtA_idxval[imin][index]), i, idxval->store[j].val * idxval->store[jj].val);
				if (debug) {
					printf("add to imin= %d index= %d i= %d val= %g\n",
					       imin, index, i, idxval->store[j].val * idxval->store[jj].val);
				}
			}
		}
		GMRFLib_idxval_free(idxval);
	}

	for (i = 0; i < g->n; i++) {
		GMRFLib_idxval_prune(AtA_idxval[i][0]);
		for (jj = 0; jj < g->lnnbs[i]; jj++) {
			GMRFLib_idxval_prune(AtA_idxval[i][1 + jj]);
		}
	}

	arg->AtA_idxval = AtA_idxval;
	arg->nlike = nlike;
	arg->like_graph = g;
	arg->like_c = Calloc(GMRFLib_MAX_THREADS, double *);
	arg->like_b = Calloc(GMRFLib_MAX_THREADS, double *);
	arg->total_const = Calloc(GMRFLib_MAX_THREADS, double);
	arg->total_b = Calloc(GMRFLib_MAX_THREADS, double *);
	for (i = 0; i < GMRFLib_MAX_THREADS; i++) {
		arg->total_b[i] = Calloc(N, double);
	}

	(*preopt)->like_graph = g;
	(*preopt)->like_Qfunc_arg = (void *) arg;
	(*preopt)->like_Qfunc = GMRFLib_preopt_like_Qfunc;
	(*preopt)->bfunc = bfunc;
	(*preopt)->like_c = arg->like_c;		       /* yes, just a copy */
	(*preopt)->like_b = arg->like_b;		       /* yes, just a copy */
	(*preopt)->total_b = arg->total_b;		       /* yes, just a copy */
	(*preopt)->total_const = arg->total_const;	       /* yes, just a copy */

	(*preopt)->n = arg->n = (*preopt)->latent_graph->n;
	(*preopt)->nlike = arg->nlike = (*preopt)->like_graph->n;
	(*preopt)->nf = arg->nf = (*preopt)->n - nbeta;
	(*preopt)->nbeta = arg->nbeta = nbeta;

	GMRFLib_graph_tp *g_arr[2];
	g_arr[0] = arg->latent_graph;
	g_arr[1] = arg->like_graph;
	GMRFLib_graph_union(&((*preopt)->preopt_graph), g_arr, 2);

	(*preopt)->preopt_Qfunc = GMRFLib_preopt_Qfunc;
	(*preopt)->preopt_Qfunc_arg = (void *) *preopt;

	GMRFLib_openmp_implement_strategy(GMRFLib_OPENMP_PLACES_DEFAULT, NULL, NULL);

	for (i = 0; i < nf; i++) {
		Free(ww[i]);
	}
	Free(ww);

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

double GMRFLib_preopt_latent_Qfunc(int node, int nnode, double *UNUSED(values), void *arg)
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

double GMRFLib_preopt_like_Qfunc(int node, int nnode, double *UNUSED(values), void *arg)
{
	if (node >= 0 && nnode < 0) {
		return NAN;
	}

	/*
	 * this is Qfunction for the likelihood part in preopt
	 */

	GMRFLib_preopt_arg_tp *a = (GMRFLib_preopt_arg_tp *) arg;

	int id = GMRFLib_thread_id, idx, k, kk, imin, imax;
	double value = 0.0, val;

	imin = IMIN(node, nnode);
	imax = IMAX(node, nnode);

	P(id);
	assert(a->like_c[id]);
	assert(a->like_b[id]);


	if (imin == imax) {
		for (kk = 0; kk < a->AtA_idxval[imin][0]->n; kk++) {
			idx = a->AtA_idxval[imin][0]->store[kk].idx;
			val = a->AtA_idxval[imin][0]->store[kk].val;
			value += val * a->like_c[id][idx];
		}
	} else {
		k = 1 + GMRFLib_iwhich_sorted(imax, a->like_graph->lnbs[imin], a->like_graph->lnnbs[imin]);
		if (k > 0) {
			for (kk = 0; kk < a->AtA_idxval[imin][k]->n; kk++) {
				idx = a->AtA_idxval[imin][k]->store[kk].idx;
				val = a->AtA_idxval[imin][k]->store[kk].val;
				value += val * a->like_c[id][idx];
			}
		} else {
			assert(k > 0);
		}
	}

	return value;
}

double GMRFLib_preopt_Qfunc(int node, int nnode, double *UNUSED(values), void *arg)
{
	if (node >= 0 && nnode < 0) {
		return NAN;
	}

	GMRFLib_preopt_tp *a = (GMRFLib_preopt_tp *) arg;
	double value = 0.0;
	int imin, imax, diag;

	imin = IMIN(node, nnode);
	imax = IMAX(node, nnode);
	diag = (imin == imax);

	if (diag || GMRFLib_graph_is_nb(imin, imax, a->like_graph)) {
		value += a->like_Qfunc(imin, imax, NULL, a->like_Qfunc_arg);
	}
	if (diag || GMRFLib_graph_is_nb(imin, imax, a->latent_graph)) {
		value += a->latent_Qfunc(imin, imax, NULL, a->latent_Qfunc_arg);
	}

	return value;
}

int GMRFLib_preopt_bnew(double *b, double *constant, GMRFLib_preopt_tp * preopt)
{
	// just add to 'b'.

	GMRFLib_preopt_bnew_latent(b, constant, preopt->n, preopt->bfunc);
	GMRFLib_preopt_bnew_like(b, preopt->like_b[GMRFLib_thread_id], preopt);

	return GMRFLib_SUCCESS;
}

int GMRFLib_preopt_bnew_latent(double *bnew, double *constant, int n, GMRFLib_bfunc_tp ** bfunc)
{
	// this is just copy of GMRFLib_bnew(...) without the alloc

	double con = 0.0, con_add = 0.0;
	if (bfunc) {
		for (int i = 0; i < n; i++) {
			if (bfunc[i]) {
				bnew[i] += GMRFLib_bfunc_eval(&con_add, bfunc[i]);
				con += con_add;
			}
		}
	}
	*constant = -con / 2.0;

	return GMRFLib_SUCCESS;
}

int GMRFLib_preopt_bnew_like(double *bnew, double *blike, GMRFLib_preopt_tp * arg)
{
	GMRFLib_preopt_arg_tp *a = (GMRFLib_preopt_arg_tp *) (arg->like_Qfunc_arg);

	for (int i = 0; i < arg->n; i++) {
		if (a->At_idxval[i]) {
			for (int jj = 0; jj < a->At_idxval[i]->n; jj++) {
				int j = a->At_idxval[i]->store[jj].idx;
				double val = a->At_idxval[i]->store[jj].val;
				bnew[i] += blike[j] * val;
			}
		}
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_preopt_update(GMRFLib_preopt_tp * preopt, double *like_b, double *like_c)
{
	int id = GMRFLib_thread_id;
	double constant = 0.0;

	printf("set like_c[%1d]\n", id);

	preopt->like_c[id] = like_c;
	preopt->like_b[id] = like_b;

	memset(preopt->total_b[id], 0, preopt->n * sizeof(double));
	GMRFLib_preopt_bnew(preopt->total_b[id], &constant, preopt);
	preopt->total_const[id] = constant;

	return GMRFLib_SUCCESS;
}

int GMRFLib_free_preopt(GMRFLib_preopt_tp * preopt)
{
	if (!preopt) {
		return GMRFLib_SUCCESS;
	}

	int i, jj;
	GMRFLib_preopt_arg_tp *a = (GMRFLib_preopt_arg_tp *) preopt->latent_Qfunc_arg;

	for (i = 0; i < preopt->n; i++) {
		GMRFLib_idxval_free(a->AtA_idxval[i][0]);
		for (jj = 0; jj < a->like_graph->lnnbs[i]; jj++) {
			GMRFLib_idxval_free(a->AtA_idxval[i][1 + jj]);
		}
	}
	Free(a->AtA_idxval);

	Free(a->idx_map_f);
	Free(a->idx_map_beta);
	Free(a->what_type);
	Free(a);

	GMRFLib_graph_free(preopt->preopt_graph);
	GMRFLib_graph_free(preopt->like_graph);
	GMRFLib_graph_free(preopt->latent_graph);
	GMRFLib_free_constr(preopt->latent_constr);
	Free(preopt);

	return GMRFLib_SUCCESS;
}

int GMRFLib_preopt_test(GMRFLib_preopt_tp * preopt)
{
	if (!preopt) {
		return GMRFLib_SUCCESS;
	}

	FIXME("LATENT");
	// GMRFLib_printf_graph(stdout, preopt->latent_graph);
	// GMRFLib_printf_Qfunc(stdout, preopt->latent_graph, preopt->latent_Qfunc, preopt->latent_Qfunc_arg);
	// GMRFLib_printf_Qfunc2(stdout, preopt->latent_graph, preopt->latent_Qfunc, preopt->latent_Qfunc_arg);
	// GMRFLib_printf_constr(stdout, preopt->latent_constr, preopt->latent_graph); 

	FIXME("LIKE");
	// GMRFLib_printf_graph(stdout, preopt->like_graph);
	// GMRFLib_printf_Qfunc(stdout, preopt->like_graph, preopt->like_Qfunc, preopt->like_Qfunc_arg);
	// GMRFLib_printf_Qfunc2(stdout, preopt->like_graph, preopt->like_Qfunc, preopt->like_Qfunc_arg);

	FIXME("JOINT");
	// GMRFLib_printf_graph(stdout, preopt->preopt_graph);
	// GMRFLib_printf_Qfunc(stdout, preopt->preopt_graph, preopt->preopt_Qfunc, preopt->preopt_Qfunc_arg);
	// GMRFLib_printf_Qfunc2(stdout, preopt->preopt_graph, preopt->preopt_Qfunc, preopt->preopt_Qfunc_arg);

	int k;
#pragma omp parallel for private(k)
	for (k = 0; k < 100; k++) {
		GMRFLib_thread_id = omp_get_thread_num();
		printf("omp id= %d\n", GMRFLib_thread_id);

		double *bb = Calloc(preopt->nlike, double);
		double *cc = Calloc(preopt->nlike, double);
		for (int i = 0; i < preopt->nlike; i++) {
			cc[i] = exp(GMRFLib_uniform());
			bb[i] = exp(GMRFLib_uniform());
		}

		GMRFLib_preopt_update(preopt, bb, cc);
		P(preopt->total_const[GMRFLib_thread_id]);
		Free(bb);
		Free(cc);
	}

	return GMRFLib_SUCCESS;
}

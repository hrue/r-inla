
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
			int npred, int nf, int **c, double **w,
			GMRFLib_graph_tp ** f_graph, GMRFLib_Qfunc_tp ** f_Qfunc,
			void **f_Qfunc_arg, char *f_sumzero, GMRFLib_constr_tp ** f_constr,
			double *f_diag,
			GMRFLib_Qfunc_tp *** ff_Qfunc, void ***ff_Qfunc_arg,
			int nbeta, double **covariate, double *prior_precision, GMRFLib_bfunc_tp ** bfunc,
			GMRFLib_ai_param_tp * UNUSED(ai_par), char *pA_fnm)
{
	if (!preopt) {
		return GMRFLib_SUCCESS;
	}
	
	GMRFLib_ENTER_ROUTINE;
	
	int i, ii, j, jj, k, kk, N, *idx_map_f = NULL, *idx_map_beta = NULL, offset, index;
	int debug = 0, nrow = 0, ncol = 0;
	double val;

	GMRFLib_idxval_tp **pA_idxval = NULL;
	GMRFLib_idxval_tp **pAA_idxval = NULL;
	GMRFLib_idxval_tp **pAAt_idxval = NULL;

	GMRFLib_constr_tp *fc = NULL;

	GMRFLib_matrix_tp * pA = NULL;

	double **ww = NULL;
	ww = Calloc(nf, double *);
	for (i = 0; i < nf; i++) {
		ww[i] = Calloc(npred, double);
		for (j = 0; j < npred; j++) {
			ww[i][j] = 1.0;
		}
	}
	if (w) {
		for (i = 0; i < nf; i++) {
			if (w[i]) {
				memcpy(ww[i], w[i], npred * sizeof(double));
			}
		}
	}

	*preopt = Calloc(1, GMRFLib_preopt_tp);
	nbeta = IMAX(0, nbeta);
	nf = IMAX(0, nf);
	(*preopt)->nf = nf;
	(*preopt)->f_Qfunc = f_Qfunc;
	(*preopt)->f_Qfunc_arg = f_Qfunc_arg;
	(*preopt)->ff_Qfunc = ff_Qfunc;
	(*preopt)->ff_Qfunc_arg = ff_Qfunc_arg;
	(*preopt)->f_graph = f_graph;
	(*preopt)->f_diag = f_diag;
	(*preopt)->nbeta = nbeta;
	(*preopt)->covariate = covariate;
	(*preopt)->prior_precision = prior_precision;

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

	(*preopt)->latent_graph = (*preopt)->latent_graph;     /* just a copy */
	(*preopt)->latent_Qfunc = GMRFLib_preopt_latent_Qfunc;
	(*preopt)->latent_Qfunc_arg = (void *) *preopt;

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

	(*preopt)->idx_map_f = idx_map_f;
	(*preopt)->idx_map_beta = idx_map_beta;
	(*preopt)->what_type = Calloc(N, GMRFLib_preopt_type_tp);
	for (i = 0; i < N; i++) {
		(*preopt)->what_type[i] = GMRFLib_preopt_what_type(i, *preopt);
	}

	if (debug) {
		printf("\tnpred %1d nf %1d nbeta %1d\n", npred, nf, nbeta);
		for (i = 0; i < npred; i++) {
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

	GMRFLib_idxval_tp **A_idxval = GMRFLib_idxval_ncreate(npred);
	for (i = 0; i < npred; i++) {
		int idx;

		for (jj = 0; jj < nf; jj++) {
			if (c[jj][i] >= 0 && ww[jj][i]) {
				idx = c[jj][i] + idx_map_f[jj];
				val = ww[jj][i];
				GMRFLib_idxval_add(&(A_idxval[i]), idx, val);
			}
		}
		for (jj = 0; jj < nbeta; jj++) {
			if (covariate[jj][i]) {
				idx = idx_map_beta[jj];
				val = covariate[jj][i];
				GMRFLib_idxval_add(&(A_idxval[i]), idx, val);
			}
		}
	}
	GMRFLib_idxval_nsort(A_idxval, npred);

	// need also At_.. below, if (pA)
	GMRFLib_idxval_tp **At_idxval = GMRFLib_idxval_ncreate(N);
	for(i = 0; i < npred; i++){
		for(k = 0; k < A_idxval[i]->n; k++){
			j = A_idxval[i]->store[k].idx;
			val = A_idxval[i]->store[k].val;
			GMRFLib_idxval_add(&(At_idxval[j]), i, val);
		}
	}
	GMRFLib_idxval_nsort(At_idxval, N);

	if (debug) {
		for(i = 0; i < npred; i++)
			GMRFLib_idxval_printf(stdout, A_idxval[i], "A_idxval");
		for(i = 0; i < N; i++)
			GMRFLib_idxval_printf(stdout, At_idxval[i], "At_idxval");
	}

	if (pA_fnm) {
		pA = GMRFLib_read_fmesher_file(pA_fnm, (long int) 0, -1);
		assert(pA);

		//if (debug) {
		//printf("read pA from [%s]\n", pA_fnm);
		//	printf("\tnrow %d ncol %d nelms %d\n", pA->nrow, pA->ncol, pA->elems);
		//for(i = 0; i < pA->elems; i++){
		//printf("\ti j x %d %d %f\n", pA->i[i], pA->j[i], pA->values[i]);
		//}
		//}

		nrow = pA->nrow;
		ncol = pA->ncol;
		assert(ncol == npred);

		// this is need to compute the linear predictor later
		pA_idxval = GMRFLib_idxval_ncreate(nrow);
		for(k = 0; k < pA->elems; k++){
			i = pA->i[k];
			j = pA->j[k];
			val = pA->values[k];
			GMRFLib_idxval_add(&(pA_idxval[i]), j, val);
		}
		GMRFLib_idxval_nsort(pA_idxval, nrow);

		GMRFLib_idx_tp **pAA_pattern = GMRFLib_idx_ncreate(nrow);
		double *row = Calloc(ncol, double);

		for(i = 0; i < nrow; i++){
			GMRFLib_matrix_get_row(row, i, pA);
			for(k = 0; k < N; k++) {
				for(jj = 0; jj < At_idxval[k]->n; jj++){
					j = At_idxval[k]->store[jj].idx;
					if (row[j]) {
						GMRFLib_idx_add(&(pAA_pattern[i]), k);
					}
				}
			}
		}
		GMRFLib_idx_nuniq(pAA_pattern, nrow);	       /* this also sorts idx's */

		if (debug) {
			char * crow = Calloc(N+1, char);
			crow[N] = '\0';
			for(i = 0; i < nrow; i++){
				memset(crow, ' ', N * sizeof(char));
				for(k = 0; k < pAA_pattern[i]->n; k++){
					j = pAA_pattern[i]->idx[k];
					//printf("Add crow i j %d %d\n", i, j);
					crow[j] = '.';
				}
				printf("pA%2d [%s]\n", i, crow);
			}
			Free(crow);
		}

		pAA_idxval = GMRFLib_idxval_ncreate(nrow);
		
		// first make a empty one filled with zeros to get the pattern. since pAA_pattern is sorted, then this will be sorted as well
		for(i = 0; i < nrow; i++){
			for(k = 0; k < pAA_pattern[i]->n; k++){
				j = pAA_pattern[i]->idx[k];
				GMRFLib_idxval_add(&(pAA_idxval[i]), j, 0.0);
			}
		}
		// then add and accumate terms using '..._addto'
		for(i = 0; i < nrow; i++){
			GMRFLib_matrix_get_row(row, i, pA);
			for(k = 0; k < N; k++) {
				for(jj = 0; jj < At_idxval[k]->n; jj++){
					j = At_idxval[k]->store[jj].idx;
					val = At_idxval[k]->store[jj].val;
					if (row[j]) {
						GMRFLib_idxval_addto(&(pAA_idxval[i]), k, row[j] * val);
					}
				}
			}
		}

		pAAt_idxval = GMRFLib_idxval_ncreate(N);
		for(i = 0; i < nrow; i++){
			for(k = 0; k < pAA_idxval[i]->n; k++){
				j = pAA_idxval[i]->store[k].idx;
				val = pAA_idxval[i]->store[k].val;
				GMRFLib_idxval_add(&(pAAt_idxval[j]), i, val);
			}
		}
		GMRFLib_idxval_nsort(pAAt_idxval, N);

		if (debug) {
			for(i = 0; i < nrow; i++){
				P(i);
				GMRFLib_idxval_printf(stdout, pAA_idxval[i], "pAA");
			}
			for(i = 0; i < N; i++){
				P(i);
				GMRFLib_idxval_printf(stdout, pAAt_idxval[i], "pAAt");
			}
		}


		for(i = 0; i < nrow; i++) {
			GMRFLib_idx_free(pAA_pattern[i]);
		}
		Free(pAA_pattern);
		Free(row);
		GMRFLib_matrix_free(pA);
	}

	// setup dimensions, see pre-opt.h for details
	if (pA_fnm) {
		(*preopt)->mpred = nrow;
		(*preopt)->npred = npred;
		(*preopt)->mnpred = npred + nrow;
		(*preopt)->Npred = nrow;
		(*preopt)->n = N;
	} else {
		(*preopt)->mpred = 0;
		(*preopt)->npred = npred;
		(*preopt)->mnpred = npred;
		(*preopt)->Npred = npred;
		(*preopt)->n = N;
	}
	
	if (debug) {
		P((*preopt)->mpred);
		P((*preopt)->npred);
		P((*preopt)->mnpred);
		P((*preopt)->Npred);
		P((*preopt)->n);
	}

	// have to create AtA from "At" & "A". the matrix 'AtA' is for the likelihood only and is either "At %*% A", or "pAAt %*% pAA",
	// depending if "pA" is there or not

	GMRFLib_idxval_tp ** gen_At;
	GMRFLib_idxval_tp ** gen_A;
	int gen_len_At;
	
	if (pA_fnm) {
		gen_A = pAA_idxval;
		gen_At = pAAt_idxval;
		gen_len_At = (*preopt)->n;
	} else {
		gen_A = A_idxval;
		gen_At = At_idxval;
		gen_len_At = (*preopt)->n;
	}

	ged = NULL;
	GMRFLib_ged_init(&ged, NULL);
	for (k = 0; k < N; k++) {
		GMRFLib_ged_add(ged, k, k);
	}

	for(i = 0; i < gen_len_At; i++){
		for (kk = 0;  kk < gen_At[i]->n; kk++) {
			k = gen_At[i]->store[kk].idx;
			for(jj = 0; jj < gen_A[k]->n; jj++){
				j = gen_A[k]->store[jj].idx;
				GMRFLib_ged_add(ged, i, j);
			}
		}
	}
	GMRFLib_graph_tp *g = NULL;
	GMRFLib_ged_build(&g, ged);
	GMRFLib_ged_free(ged);
	assert(g->n == gen_len_At);

	GMRFLib_idxval_tp ***AtA_idxval = Calloc(gen_len_At, GMRFLib_idxval_tp **);
	for (i = 0; i < g->n; i++) {
		AtA_idxval[i] = GMRFLib_idxval_ncreate(1 + g->lnnbs[i]);
	}

	for(i = 0; i < gen_len_At; i++) {
		for (kk = 0;  kk < gen_At[i]->n; kk++) {
			k = gen_At[i]->store[kk].idx;
			for(jj = 0; jj < gen_A[k]->n; jj++){
				j = gen_A[k]->store[jj].idx;
				if (j >= i) {
					if (i == j) {
						index = 0;
					} else {
						index = 1 + GMRFLib_iwhich_sorted(j, g->lnbs[i], g->lnnbs[i]);
						assert(index > 0);
					}
					GMRFLib_idxval_add(&(AtA_idxval[i][index]), k,
							   gen_At[i]->store[kk].val * gen_A[k]->store[jj].val);
				}
			}
		}
	}

	if (debug) {
		FIXME("AtA");
		for(i = 0; i < N; i++) {
			double sum = 0.0;
			
			printf("term %d %d\n", i, i);
			for(kk = 0; kk < AtA_idxval[i][0]->n; kk++){
				printf("\tkk idx val %d %d %f\n", kk,
					AtA_idxval[i][0]->store[kk].idx,
					AtA_idxval[i][0]->store[kk].val);
				sum += AtA_idxval[i][0]->store[kk].val;
			}
			printf("\tsum %g\n", sum);

			for(jj= 0; jj < g->lnnbs[i]; jj++){
				j = g->lnbs[i][jj];
				printf("term %d %d\n", i, j);
				sum = 0.0;
				for(kk = 0; kk < AtA_idxval[i][1+jj]->n; kk++){
					printf("\tkk idx val %d %d %f\n", kk,
						AtA_idxval[i][1+jj]->store[kk].idx,
						AtA_idxval[i][1+jj]->store[kk].val);
					sum += AtA_idxval[i][1+jj]->store[kk].val;
				}
				printf("\tsum %g\n", sum);
			}
		}
	}
	
	GMRFLib_idxval_nprune(pA_idxval, nrow);
	GMRFLib_idxval_nprune(pAAt_idxval, N);
	GMRFLib_idxval_nprune(A_idxval, npred);
	GMRFLib_idxval_nprune(At_idxval, N);
	for (i = 0; i < g->n; i++) {
		GMRFLib_idxval_nprune(AtA_idxval[i], 1 + g->lnnbs[i]);
	}

	(*preopt)->A_idxval = A_idxval;
	(*preopt)->At_idxval = At_idxval;
	(*preopt)->pA_idxval = pA_idxval;
	(*preopt)->pAAt_idxval = pAAt_idxval;
	(*preopt)->AtA_idxval = AtA_idxval;

	(*preopt)->like_graph = g;
	(*preopt)->like_c = Calloc(GMRFLib_MAX_THREADS, double *);
	(*preopt)->like_b = Calloc(GMRFLib_MAX_THREADS, double *);
	(*preopt)->total_b = Calloc(GMRFLib_MAX_THREADS, double *);
	(*preopt)->total_const = Calloc(GMRFLib_MAX_THREADS, double);

	(*preopt)->like_Qfunc_arg = (void *) *preopt;
	(*preopt)->like_Qfunc = GMRFLib_preopt_like_Qfunc;
	(*preopt)->bfunc = bfunc;
	(*preopt)->nf = (*preopt)->n - nbeta;
	(*preopt)->nbeta = nbeta;

	GMRFLib_graph_tp *g_arr[2];
	g_arr[0] = (*preopt)->latent_graph;
	g_arr[1] = (*preopt)->like_graph;
	GMRFLib_graph_union(&((*preopt)->preopt_graph), g_arr, 2);

	(*preopt)->preopt_Qfunc = GMRFLib_preopt_Qfunc;
	(*preopt)->preopt_Qfunc_arg = (void *) *preopt;

	GMRFLib_openmp_implement_strategy(GMRFLib_OPENMP_PLACES_DEFAULT, NULL, NULL);

	if (pAA_idxval) {
		for(i = 0; i < nrow; i++) {
			GMRFLib_idxval_free(pAA_idxval[i]);
		}
		Free(pAA_idxval);
	}

	for (i = 0; i < nf; i++) {
		Free(ww[i]);
	}
	Free(ww);

	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}

GMRFLib_preopt_type_tp GMRFLib_preopt_what_type(int node, GMRFLib_preopt_tp * preopt)
{
	int i;
	GMRFLib_preopt_type_tp t = { GMRFLib_PREOPT_TP___VOID, -1, -1 };

	if (preopt->nf && node < preopt->idx_map_f[preopt->nf]) {
		t.tp = GMRFLib_PREOPT_TP_F;
		for (i = 0; i < preopt->nf; i++) {
			if (node < preopt->idx_map_f[i + 1]) {
				t.tp_idx = i;
				t.idx = node - preopt->idx_map_f[i];
				break;
			}
		}
	} else if (preopt->nbeta && node < preopt->idx_map_beta[preopt->nbeta]) {
		t.tp = GMRFLib_PREOPT_TP_BETA;
		t.tp_idx = node - preopt->idx_map_beta[0];
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
	GMRFLib_preopt_tp *a = NULL;
	int ii, jj;
	double value = 0.0;
	GMRFLib_preopt_type_tp it, jt;

	a = (GMRFLib_preopt_tp *) arg;
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
				if (it.idx == jt.idx) {
					value += a->f_diag[it.tp_idx];
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

	GMRFLib_preopt_tp *a = (GMRFLib_preopt_tp *) arg;

	int id = GMRFLib_thread_id, idx, k, kk, imin, imax;
	double value = 0.0, val;

	imin = IMIN(node, nnode);
	imax = IMAX(node, nnode);

	if (a->like_c[id]) {
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
	} else {
		value = 0.0;
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

	// GMRFLib_preopt_bnew_latent(b, constant, preopt->n, preopt->bfunc);
	*constant = 0.0;

	GMRFLib_preopt_bnew_like(b, preopt->like_b[GMRFLib_thread_id], preopt);

	return GMRFLib_SUCCESS;
}

int GMRFLib_preopt_bnew_latent(double *bnew, double *constant, int n, GMRFLib_bfunc_tp ** bfunc)
{
	FIXME("NOT IN USE");
	abort();

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

int GMRFLib_preopt_bnew_like(double *bnew, double *blike, GMRFLib_preopt_tp * preopt)
{
	// add to 'bnew' the contribution from the likelihood for preopt

	GMRFLib_idxval_tp **A = NULL;
	if (preopt->pA_idxval) {
		A = preopt->pAAt_idxval;
		assert(preopt->mpred > 0);
	} else {
		A = preopt->At_idxval;
  		assert(preopt->mpred == 0);
	}

	for (int i = 0; i < preopt->n; i++) {
		double val;
		int idx;
		if (A[i]) {
			for (int jj = 0; jj < A[i]->n; jj++) {
				idx = A[i]->store[jj].idx;
				val = A[i]->store[jj].val;
				bnew[i] += blike[idx] * val;
			}
		}
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_preopt_predictor(double *predictor, double *latent, GMRFLib_preopt_tp * preopt) 
{
	return GMRFLib_preopt_predictor_core(predictor, latent, preopt, 1);
}

int GMRFLib_preopt_full_predictor(double *predictor, double *latent, GMRFLib_preopt_tp * preopt)
{
	return GMRFLib_preopt_predictor_core(predictor, latent, preopt, 0);
}

int GMRFLib_preopt_predictor_core(double *predictor, double *latent, GMRFLib_preopt_tp * preopt,
				  int likelihood_only)
{
	// if likelihood_only, only compute the part that is needed for the likelihood.

	// if !likelihood_only, compute the whole likelihood

	int offset = 0;
	double *pred = Calloc(preopt->mnpred, double);

	if (preopt->pA_idxval) {
		offset = preopt->mpred;
		assert(preopt->mpred > 0);
	} else {
  		assert(preopt->mpred == 0);
	}

	for (int i = 0; i < preopt->npred; i++) {
		if (preopt->A_idxval[i]) {
			for (int jj = 0; jj < preopt->A_idxval[i]->n; jj++) {
				int idx = preopt->A_idxval[i]->store[jj].idx;
				double val = preopt->A_idxval[i]->store[jj].val;
				pred[offset + i] += val * latent[idx];
			}
		}
	}
	for (int i = 0; i < preopt->mpred; i++) {
		if (preopt->pA_idxval[i]) {
			for (int jj = 0; jj < preopt->pA_idxval[i]->n; jj++) {
				int idx = preopt->pA_idxval[i]->store[jj].idx;
				double val = preopt->pA_idxval[i]->store[jj].val;
				pred[i] += val * pred[offset + idx];
			}
		}
	}

	if (likelihood_only) {
		memcpy(predictor, pred, preopt->Npred * sizeof(double));
	} else {
		memcpy(predictor, pred, preopt->mnpred * sizeof(double));
	}
	Free(pred);

	return GMRFLib_SUCCESS;
}

int GMRFLib_preopt_update(GMRFLib_preopt_tp * preopt, double *like_b, double *like_c)
{
	int id = GMRFLib_thread_id, np = preopt->Npred;
	double constant = 0.0;

	if (!(preopt->like_b[id])) {
		preopt->like_b[id] = Calloc(np, double);
		preopt->like_c[id] = Calloc(np, double);
	}
	memcpy(preopt->like_b[id], like_b, np * sizeof(double));
	memcpy(preopt->like_c[id], like_c, np * sizeof(double));

	if (!(preopt->total_b[id])) {
		preopt->total_b[id] = Calloc(preopt->n, double);
	}
	memset(preopt->total_b[id], 0, preopt->n * sizeof(double));

	GMRFLib_preopt_bnew(preopt->total_b[id], &constant, preopt);
	preopt->total_const[id] = constant;

	return GMRFLib_SUCCESS;
}

int GMRFLib_preopt_free(GMRFLib_preopt_tp * preopt)
{
	if (!preopt) {
		return GMRFLib_SUCCESS;
	}

	int i, jj;

	for (i = 0; i < preopt->n; i++) {
		GMRFLib_idxval_free(preopt->AtA_idxval[i][0]);
		for (jj = 0; jj < preopt->like_graph->lnnbs[i]; jj++) {
			GMRFLib_idxval_free(preopt->AtA_idxval[i][1 + jj]);
		}
	}
	Free(preopt->AtA_idxval);

	if (preopt->pA_idxval) {
		for(i = 0; i < preopt->mpred; i++) {
			GMRFLib_idxval_free(preopt->pA_idxval[i]);
		}
		for(i = 0; i < preopt->n; i++) {
			GMRFLib_idxval_free(preopt->pAAt_idxval[i]);
		}
	}
	for(i = 0; i < preopt->npred; i++) {
		GMRFLib_idxval_free(preopt->A_idxval[i]);
	}
	for(i = 0; i < preopt->n; i++) {
		GMRFLib_idxval_free(preopt->At_idxval[i]);
	}

	Free(preopt->idx_map_f);
	Free(preopt->idx_map_beta);
	Free(preopt->what_type);

	for (i = 0; i < GMRFLib_MAX_THREADS; i++) {
		Free(preopt->like_b[i]);
		Free(preopt->like_c[i]);
		Free(preopt->total_b[i]);
	}
	Free(preopt->like_b);
	Free(preopt->like_c);
	Free(preopt->total_b);

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
	GMRFLib_printf_Qfunc2(stdout, preopt->latent_graph, preopt->latent_Qfunc, preopt->latent_Qfunc_arg);
	// GMRFLib_printf_constr(stdout, preopt->latent_constr, preopt->latent_graph); 

	FIXME("LIKE");
	// GMRFLib_printf_graph(stdout, preopt->like_graph);
	// GMRFLib_printf_Qfunc(stdout, preopt->like_graph, preopt->like_Qfunc, preopt->like_Qfunc_arg);
	GMRFLib_printf_Qfunc2(stdout, preopt->like_graph, preopt->like_Qfunc, preopt->like_Qfunc_arg);

	FIXME("JOINT");
	// GMRFLib_printf_graph(stdout, preopt->preopt_graph);
	// GMRFLib_printf_Qfunc(stdout, preopt->preopt_graph, preopt->preopt_Qfunc, preopt->preopt_Qfunc_arg);
	GMRFLib_printf_Qfunc2(stdout, preopt->preopt_graph, preopt->preopt_Qfunc, preopt->preopt_Qfunc_arg);

	return GMRFLib_SUCCESS;
}

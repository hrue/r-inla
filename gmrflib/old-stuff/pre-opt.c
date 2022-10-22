
/* pre-opt.c
 * 
 * Copyright (C) 2021-2022 Havard Rue
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

int integer_one = 1;
#define DDOT(N_, X_, Y_) ddot_(&(N_), X_, &integer_one, Y_, &integer_one)

#define DOT_PRODUCT_GROUP(VALUE_, ELM_, ARR_)				\
	if (1) {							\
		double value_ = 0.0;					\
		for (int g_ = 0; g_ < ELM_->g_n; g_++) {		\
			int istart_ = ELM_->g_i[g_];			\
			int *ii_ = &(ELM_->idx[istart_]);		\
			double *vv_ = &(ELM_->val[istart_]);		\
			switch(ELM_->g_len[g_])				\
			{						\
			case 1:						\
				value_ +=				\
					vv_[0] * ARR_[ii_[0]];		\
				break;					\
			case 2:						\
				value_ +=				\
					vv_[0] * ARR_[ii_[0]] +		\
					vv_[1] * ARR_[ii_[1]];		\
				break;					\
			case 3:						\
				value_ +=				\
					vv_[0] * ARR_[ii_[0]] +		\
					vv_[1] * ARR_[ii_[1]] +		\
					vv_[2] * ARR_[ii_[2]];		\
				break;					\
			default:					\
				value_ += DDOT(ELM_->g_len[g_], vv_, &(ARR_[ii_[0]])); \
			}						\
		}							\
		VALUE_ = value_;					\
	}

#define DOT_PRODUCT_SERIAL(VALUE_, ELM_, ARR_)				\
	if (1) {							\
		double value_ = 0.0;					\
		_Pragma("GCC ivdep")					\
			_Pragma("GCC unroll 8")				\
			for (int i_= 0; i_ < ELM_->n; i_++) {		\
				value_ += ELM_->val[i_] * ARR_[ELM_->idx[i_]]; \
			}						\
		VALUE_ = value_;					\
	}

#define DOT_PRODUCT(VALUE_, ELM_, ARR_)			\
	if (GMRFLib_preopt_like_strategy == 0) {	\
		DOT_PRODUCT_SERIAL(VALUE_, ELM_, ARR_);	\
	} else {					\
		DOT_PRODUCT_GROUP(VALUE_, ELM_, ARR_);	\
	}

int GMRFLib_preopt_init(GMRFLib_preopt_tp ** preopt,
			int npred, int nf, int **c, double **w,
			GMRFLib_graph_tp ** f_graph, GMRFLib_Qfunc_tp ** f_Qfunc,
			void **f_Qfunc_arg, char *f_sumzero, GMRFLib_constr_tp ** f_constr,
			double *f_diag,
			GMRFLib_Qfunc_tp *** ff_Qfunc, void ***ff_Qfunc_arg,
			int nbeta, double **covariate, double *prior_precision, GMRFLib_bfunc_tp ** bfunc,
			GMRFLib_ai_param_tp * UNUSED(ai_par), char *pA_fnm)
{
#define SHOW_TIME(_msg)							\
	if (debug) {							\
		printf("\t\tGMRFLib_preopt_init: %-16s %7.2fs\n", _msg, GMRFLib_cpu() - tref); \
		tref =  GMRFLib_cpu();					\
	}

	if (!preopt) {
		return GMRFLib_SUCCESS;
	}
	GMRFLib_ENTER_ROUTINE;

	int i, ii, j, jj, k, kk, N = 0, *idx_map_f = NULL, *idx_map_beta = NULL, offset, index;
	int nrow = 0, ncol = 0;
	int debug = GMRFLib_DEBUG_IF_TRUE();
	int debug_detailed = 0;
	int id = GMRFLib_thread_id;

	double tref = GMRFLib_cpu();
	double **ww = NULL;
	GMRFLib_constr_tp *fc = NULL;
	GMRFLib_idx_tp **pAA_pattern = NULL;
	GMRFLib_idxval_tp **A_idxval = NULL;
	GMRFLib_idxval_tp ***AtA_idxval = NULL;
	GMRFLib_idxval_tp **At_idxval = NULL;
	GMRFLib_idxval_tp **pAA_idxval = NULL;
	GMRFLib_idxval_tp **pAAt_idxval = NULL;
	GMRFLib_idxval_tp **pA_idxval = NULL;
	GMRFLib_matrix_tp *pA = NULL;
	GMRFLib_idxval_tp *elm = NULL;

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
				Memcpy(ww[i], w[i], npred * sizeof(double));
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

	SHOW_TIME("setup1");

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
	N = offset;					       /* N is the size of the latent (no predictors) */

	/*
	 * If we have cross-terms, make sure to mark these cross-terms as neigbours; just fill them with zero's.
	 * If they are connected in the data, then this value
	 * will be overrided and this is how it should be.
	 */
	GMRFLib_ged_tp *ged = NULL;
	GMRFLib_ged_init2(&ged, N);

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
	(*preopt)->latent_Qfunc = GMRFLib_preopt_latent_Qfunc;
	(*preopt)->latent_Qfunc_arg = (void *) *preopt;

	/*
	 * build the constraint, if any. 
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

	if (debug_detailed) {
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

	SHOW_TIME("setup2");

	// build up structure for the likelihood part

	int nt = -1;
	if (GMRFLib_OPENMP_IN_PARALLEL()) {
		nt = GMRFLib_openmp->max_threads_inner;
	} else {
		nt = GMRFLib_openmp->max_threads_outer;
	}
	nt = IMIN(GMRFLib_MAX_THREADS(), nt);

	A_idxval = GMRFLib_idxval_ncreate(npred);
#pragma omp parallel for num_threads(nt)
	for (int i = 0; i < npred; i++) {
		GMRFLib_thread_id = id;

		int idx;
		double val;

		for (int jj = 0; jj < nf; jj++) {
			if (c[jj][i] >= 0 && ww[jj][i]) {
				idx = c[jj][i] + idx_map_f[jj];
				val = ww[jj][i];
				GMRFLib_idxval_add(&(A_idxval[i]), idx, val);
			}
		}
		for (int jj = 0; jj < nbeta; jj++) {
			if (covariate[jj][i]) {
				idx = idx_map_beta[jj];
				val = covariate[jj][i];
				GMRFLib_idxval_add(&(A_idxval[i]), idx, val);
			}
		}
		GMRFLib_idxval_sort(A_idxval[i]);
	}
	GMRFLib_idxval_to_matrix(&((*preopt)->A), A_idxval, npred, N);
	SHOW_TIME("A_idxval");


	// need also At_.. below, if (pA)
	At_idxval = GMRFLib_idxval_ncreate(N);
	for (i = 0; i < npred; i++) {
		elm = A_idxval[i];
		for (k = 0; k < A_idxval[i]->n; k++) {
			GMRFLib_idxval_add(&(At_idxval[elm->idx[k]]), i, elm->val[k]);
		}
	}
	GMRFLib_idxval_nsort(At_idxval, N, nt);

	SHOW_TIME("At_idxval");
	if (debug_detailed) {
		for (i = 0; i < npred; i++) {
			GMRFLib_idxval_printf(stdout, A_idxval[i], "A_idxval");
		}
		for (i = 0; i < N; i++) {
			GMRFLib_idxval_printf(stdout, At_idxval[i], "At_idxval");
		}
	}

	if (pA_fnm) {
		pA = GMRFLib_read_fmesher_file(pA_fnm, (long int) 0, -1);
		assert(pA);

		if (debug_detailed) {
			printf("read pA from [%s]\n", pA_fnm);
			printf("\tnrow %d ncol %d nelms %d\n", pA->nrow, pA->ncol, pA->elems);
			for (i = 0; i < pA->elems; i++) {
				printf("\ti j x %d %d %f\n", pA->i[i], pA->j[i], pA->values[i]);
			}
		}

		nrow = pA->nrow;
		ncol = pA->ncol;
		assert(ncol == npred);
		SHOW_TIME("read pA");

		// this is need to compute the linear predictor later
		pA_idxval = GMRFLib_idxval_ncreate(nrow);
		for (k = 0; k < pA->elems; k++) {
			i = pA->i[k];
			j = pA->j[k];
			GMRFLib_idxval_add(&(pA_idxval[i]), j, pA->values[k]);
		}
		GMRFLib_idxval_nsort(pA_idxval, nrow, nt);
		(*preopt)->pA = pA;
		SHOW_TIME("create pA_idxval");

		pAA_pattern = GMRFLib_idx_ncreate(nrow);
#pragma omp parallel for private (i, k, kk, j, jj) num_threads(nt)
		for (i = 0; i < nrow; i++) {
			GMRFLib_thread_id = id;
			GMRFLib_idxval_tp *row_idxval = NULL;
			GMRFLib_matrix_get_row_idxval(&row_idxval, i, pA);

			for (jj = 0; jj < row_idxval->n; jj++) {
				j = row_idxval->idx[jj];
				// find possible k's
				for (kk = 0; kk < A_idxval[j]->n; kk++) {
					k = A_idxval[j]->idx[kk];
					GMRFLib_idx_add(&(pAA_pattern[i]), k);
				}
			}
			GMRFLib_idx_uniq(pAA_pattern[i]);      /* also sorts */
			GMRFLib_idxval_free(row_idxval);
		}
		SHOW_TIME("pAA_pattern");

		if (debug_detailed) {
			char *crow = Calloc(N + 1, char);
			crow[N] = '\0';
			for (i = 0; i < nrow; i++) {
				Memset(crow, ' ', N * sizeof(char));
				for (k = 0; k < pAA_pattern[i]->n; k++) {
					j = pAA_pattern[i]->idx[k];
					// printf("Add crow i j %d %d\n", i, j);
					crow[j] = '.';
				}
				printf("pAA%2d [%s]\n", i, crow);
			}
			Free(crow);
		}
		// first make a empty one filled with zeros to get the pattern. since pAA_pattern is sorted, then this will be sorted as well
		pAA_idxval = GMRFLib_idxval_ncreate(nrow);
#pragma omp parallel for private (i, k, j) num_threads(nt)
		for (i = 0; i < nrow; i++) {
			GMRFLib_thread_id = id;
			int *idx = pAA_pattern[i]->idx;
			for (k = 0; k < pAA_pattern[i]->n; k++) {
				GMRFLib_idxval_add(&(pAA_idxval[i]), idx[k], 0.0);
			}
		}

		SHOW_TIME("init pAA_idxval");

#pragma omp parallel for private (i, k, kk, j, jj) num_threads(nt)
		for (i = 0; i < nrow; i++) {
			GMRFLib_thread_id = id;

			int step, s, ia;
			int steps[] = { 262144, 32768, 4096, 512, 64, 8, 1 };
			int nsteps = sizeof(steps) / sizeof(int);

			int row_n;
			GMRFLib_idxval_tp *row_idxval = NULL;
			GMRFLib_idxval_tp *row_elm = NULL;

			GMRFLib_matrix_get_row_idxval(&row_idxval, i, pA);
			row_elm = row_idxval;
			row_n = row_idxval->n;

			for (jj = 0; jj < pAA_pattern[i]->n; jj++) {
				j = pAA_pattern[i]->idx[jj];
				GMRFLib_idxval_tp *At_elm = At_idxval[j];

				int At_n = At_idxval[j]->n, irow = 0, iAt = 0;
				while (irow < row_n && iAt < At_n) {

					k = row_elm->idx[irow];
					kk = At_elm->idx[iAt];

					if (k < kk) {
						irow++;
						for (s = 0; s < nsteps; s++) {
							step = steps[s];
							if (step < row_n) {
								ia = irow + step;
								while (ia < row_n && row_elm->idx[ia] < kk)
									ia += step;
								irow = ia - step;
							}
						}
					} else if (k > kk) {
						iAt++;
						for (s = 0; s < nsteps; s++) {
							step = steps[s];
							if (step < At_n) {
								ia = iAt + step;
								while (ia < At_n && At_elm->idx[ia] < k)
									ia += step;
								iAt = ia - step;
							}
						}
					} else {
						GMRFLib_idxval_addto(&(pAA_idxval[i]), j, row_elm->val[irow] * At_elm->val[iAt]);
						irow++;
						iAt++;
					}
				}
			}
			GMRFLib_idxval_free(row_idxval);
		}
		SHOW_TIME("pAA_idxval");

		pAAt_idxval = GMRFLib_idxval_ncreate(N);
		for (i = 0; i < nrow; i++) {
			elm = pAA_idxval[i];
			for (k = 0; k < pAA_idxval[i]->n; k++) {
				GMRFLib_idxval_add(&(pAAt_idxval[elm->idx[k]]), i, elm->val[k]);
			}
		}
		GMRFLib_idxval_nsort(pAAt_idxval, N, 0);       /* as N is typical small */
		SHOW_TIME("pAAt_idxval");

		if (debug_detailed) {
			for (i = 0; i < nrow; i++) {
				P(i);
				GMRFLib_idxval_printf(stdout, pAA_idxval[i], "pAA");
			}
			for (i = 0; i < N; i++) {
				P(i);
				GMRFLib_idxval_printf(stdout, pAAt_idxval[i], "pAAt");
			}
		}

		for (i = 0; i < nrow; i++) {
			GMRFLib_idx_free(pAA_pattern[i]);
		}
		Free(pAA_pattern);
		SHOW_TIME("End pA... ");
	}
	// setup dimensions, see pre-opt.h for the details
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
	// we have to create AtA from "At" & "A". the matrix 'AtA' is for the likelihood only and is either "At %*% A", or "pAAt %*% pAA",
	// depending if "pA" is there or not

	GMRFLib_idxval_tp **gen_At = NULL;
	GMRFLib_idxval_tp **gen_A = NULL;
	int gen_len_At = -1;

	if (pA_fnm) {
		gen_A = pAA_idxval;
		gen_At = pAAt_idxval;
		gen_len_At = (*preopt)->n;
	} else {
		gen_A = A_idxval;
		gen_At = At_idxval;
		gen_len_At = (*preopt)->n;
	}

	SHOW_TIME("admin1");

	GMRFLib_graph_tp *g = NULL;
	ged = NULL;
	GMRFLib_ged_init2(&ged, N);
	for (i = 0; i < gen_len_At; i++) {
		for (kk = 0; kk < gen_At[i]->n; kk++) {
			k = gen_At[i]->idx[kk];
			for (jj = 0; jj < gen_A[k]->n; jj++) {
				j = gen_A[k]->idx[jj];
				GMRFLib_ged_add(ged, i, j);
			}
		}
	}
	SHOW_TIME("build graph part 1");

	GMRFLib_ged_build(&g, ged);
	GMRFLib_ged_free(ged);
	assert(g->n == gen_len_At);
	SHOW_TIME("build graph part 2");

	AtA_idxval = Calloc(gen_len_At, GMRFLib_idxval_tp **);
	for (i = 0; i < g->n; i++) {
		AtA_idxval[i] = GMRFLib_idxval_ncreate(1 + g->lnnbs[i]);
	}

#pragma omp parallel for private (i, kk, k, jj, j, index) num_threads(nt)
	for (i = 0; i < gen_len_At; i++) {
		GMRFLib_thread_id = id;
		int guess[] = { 0, 0 };
		for (kk = 0; kk < gen_At[i]->n; kk++) {
			k = gen_At[i]->idx[kk];
			for (jj = 0; jj < gen_A[k]->n; jj++) {
				j = gen_A[k]->idx[jj];
				if (j >= i) {
					if (i == j) {
						index = 0;
					} else {
						index = 1 + GMRFLib_iwhich_sorted(j, g->lnbs[i], g->lnnbs[i], guess);
						assert(index > 0);
					}
					GMRFLib_idxval_add(&(AtA_idxval[i][index]), k, gen_At[i]->val[kk] * gen_A[k]->val[jj]);
				}
			}
		}
	}
	SHOW_TIME("AtA_idxval");

	if (debug_detailed) {
		FIXME("AtA");
		for (i = 0; i < N; i++) {
			double sum = 0.0;
			printf("term %d %d\n", i, i);
			for (kk = 0; kk < AtA_idxval[i][0]->n; kk++) {
				printf("\tkk idx val %d %d %f\n", kk, AtA_idxval[i][0]->idx[kk], AtA_idxval[i][0]->val[kk]);
				sum += AtA_idxval[i][0]->val[kk];
			}
			printf("\tsum %g\n", sum);
			for (jj = 0; jj < g->lnnbs[i]; jj++) {
				j = g->lnbs[i][jj];
				printf("term %d %d\n", i, j);
				sum = 0.0;
				for (kk = 0; kk < AtA_idxval[i][1 + jj]->n; kk++) {
					printf("\tkk idx val %d %d %f\n", kk, AtA_idxval[i][1 + jj]->idx[kk], AtA_idxval[i][1 + jj]->val[kk]);
					sum += AtA_idxval[i][1 + jj]->val[kk];
				}
				printf("\tsum %g\n", sum);
			}
		}
	}
#pragma omp parallel for private (i) num_threads(nt)
	for (i = 0; i < g->n; i++) {
		GMRFLib_thread_id = id;
		GMRFLib_idxval_nsort(AtA_idxval[i], 1 + g->lnnbs[i], 0);
	}
	SHOW_TIME("sort AtA_idxval");

	(*preopt)->A_idxval = A_idxval;
	(*preopt)->At_idxval = At_idxval;
	(*preopt)->pA_idxval = pA_idxval;
	(*preopt)->pAA_idxval = pAA_idxval;
	(*preopt)->pAAt_idxval = pAAt_idxval;
	(*preopt)->AtA_idxval = AtA_idxval;

	(*preopt)->like_graph = g;
	(*preopt)->like_c = Calloc(GMRFLib_MAX_THREADS(), double *);
	(*preopt)->like_b = Calloc(GMRFLib_MAX_THREADS(), double *);
	(*preopt)->total_b = Calloc(GMRFLib_MAX_THREADS(), double *);


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

	for (i = 0; i < nf; i++) {
		Free(ww[i]);
	}
	Free(ww);

	GMRFLib_openmp_implement_strategy(GMRFLib_OPENMP_PLACES_DEFAULT, NULL, NULL);

	SHOW_TIME("admin2");
#undef  SHOW_TIME

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
	// as this one is always called through preopt_Qfunc
	// assert(nnode >= node);
	// if (node >= 0 && nnode < 0) return NAN;

	/*
	 * this is Qfunction for the preopt-function 
	 */
	GMRFLib_preopt_tp *a = NULL;
	GMRFLib_preopt_type_tp it, jt;
	double value = 0.0;
	int same_tp, same_idx;

	a = (GMRFLib_preopt_tp *) arg;
	it = a->what_type[node];
	jt = a->what_type[nnode];
	same_tp = (it.tp_idx == jt.tp_idx);
	same_idx = (it.idx == jt.idx);

	switch (it.tp) {
	case GMRFLib_PREOPT_TP_F:
		switch (jt.tp) {
		case GMRFLib_PREOPT_TP_F:
			if (same_tp) {
				if (same_idx || GMRFLib_graph_is_nb(it.idx, jt.idx, a->f_graph[it.tp_idx])) {
					value += a->f_Qfunc[it.tp_idx] (it.idx, jt.idx, NULL, (a->f_Qfunc_arg ? a->f_Qfunc_arg[it.tp_idx] : NULL));
				}
				if (same_idx) {
					value += a->f_diag[it.tp_idx];
				}
			}

			/*
			 * only for the same index and different types; used to define `interaction between fields'. this is a 'workaround' for a INLA problem.. 
			 */
			if (a->ff_Qfunc) {
				if (same_idx && !same_tp && a->ff_Qfunc[it.tp_idx][jt.tp_idx]) {
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
			if (same_tp) {
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
	/*
	 * this is Qfunction for the likelihood part in preopt
	 */
	GMRFLib_preopt_tp *a = (GMRFLib_preopt_tp *) arg;
	GMRFLib_idxval_tp *elm = NULL;
	int id = GMRFLib_thread_id;
	double *lc = a->like_c[id], value = 0.0;

	if (!lc) {
		return 0.0;
	}
	// imin = node; imax = nnode;
	if (node == nnode) {
		elm = a->AtA_idxval[node][0];
		DOT_PRODUCT(value, elm, lc);
	} else {
		// use also this [low, high] guess, which is updated automatically
		static int guess[] = { 0, 0 };
#pragma omp threadprivate(guess)

		int k = 1 + GMRFLib_iwhich_sorted(nnode, a->like_graph->lnbs[node], a->like_graph->lnnbs[node], guess);
		elm = a->AtA_idxval[node][k];
		DOT_PRODUCT(value, elm, lc);
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

double GMRFLib_preopt_Qfunc_like(int node, int nnode, double *UNUSED(values), void *arg)
{
	// standalone function to return the likelihood part only
	if (node >= 0 && nnode < 0) {
		return NAN;
	}

	GMRFLib_preopt_tp *a = (GMRFLib_preopt_tp *) arg;
	int imin, imax, diag;

	imin = IMIN(node, nnode);
	imax = IMAX(node, nnode);
	diag = (imin == imax);

	double value = 0.0;
	if (diag || GMRFLib_graph_is_nb(imin, imax, a->like_graph)) {
		value += a->like_Qfunc(imin, imax, NULL, a->like_Qfunc_arg);
	}

	return value;
}

double GMRFLib_preopt_Qfunc_prior(int node, int nnode, double *UNUSED(values), void *arg)
{
	// standalone function to return the prior part
	if (node >= 0 && nnode < 0) {
		return NAN;
	}

	GMRFLib_preopt_tp *a = (GMRFLib_preopt_tp *) arg;
	int imin, imax, diag;

	imin = IMIN(node, nnode);
	imax = IMAX(node, nnode);
	diag = (imin == imax);

	double value = 0.0;
	if (diag || GMRFLib_graph_is_nb(imin, imax, a->latent_graph)) {
		value = a->latent_Qfunc(imin, imax, NULL, a->latent_Qfunc_arg);
	}

	return value;
}


int GMRFLib_preopt_bnew(double *b, GMRFLib_preopt_tp * preopt)
{
	GMRFLib_ENTER_ROUTINE;
	GMRFLib_preopt_bnew_like(b, preopt->like_b[GMRFLib_thread_id], preopt);

	GMRFLib_LEAVE_ROUTINE;
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

#define CODE_BLOCK						\
	for (int i = 0; i < preopt->n; i++) {			\
		CODE_BLOCK_SET_THREAD_ID();			\
		if (A[i]) {					\
			GMRFLib_idxval_tp *elm = A[i];		\
			DOT_PRODUCT(bnew[i], elm, blike);	\
		}						\
	}

	RUN_CODE_BLOCK(GMRFLib_MAX_THREADS(), 0, 0);
#undef CODE_BLOCK

	return GMRFLib_SUCCESS;
}

int GMRFLib_preopt_bnew_like_ORIG(double *bnew, double *blike, GMRFLib_preopt_tp * preopt)
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

#define CODE_BLOCK							\
	for (int i = 0; i < preopt->n; i++) {				\
		CODE_BLOCK_SET_THREAD_ID();				\
		if (A[i]) {						\
			GMRFLib_idxval_tp *elm = A[i];			\
			for (int jj = 0; jj < A[i]->n; jj++) {		\
				bnew[i] += blike[elm->idx[jj]] * elm->val[jj]; \
			}						\
		}							\
	}

	RUN_CODE_BLOCK(GMRFLib_MAX_THREADS(), 0, 0);
#undef CODE_BLOCK

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

int GMRFLib_preopt_predictor_core(double *predictor, double *latent, GMRFLib_preopt_tp * preopt, int likelihood_only)
{
	// if likelihood_only, only compute the part that is needed for the likelihood.

	// if !likelihood_only, compute the whole predictor

	GMRFLib_ENTER_ROUTINE;

	double *pred = Calloc(preopt->mnpred, double);
	int data_rich_case = (IMAX(preopt->mpred, preopt->npred) > preopt->n);
	int offset = 0;

	if (preopt->pA_idxval) {
		offset = preopt->mpred;
		assert(preopt->mpred > 0);
	} else {
		assert(preopt->mpred == 0);
	}

	if (data_rich_case) {

		// data-rich case

		if (preopt->pA_idxval) {
			// both loops
#define CODE_BLOCK							\
			for(int j = 0; j < 2; j++) {			\
				CODE_BLOCK_SET_THREAD_ID();		\
				if (j == 0) {				\
					for(int i = 0; i < preopt->n; i++) { \
						GMRFLib_idxval_tp *At = preopt->At_idxval[i]; \
						double lat = latent[i];	\
						for(int k = 0; k < At->n; k++) { \
							pred[offset + At->idx[k]] += lat * At->val[k]; \
						}			\
					}				\
				} else {				\
					for(int i = 0; i < preopt->n; i++) { \
						GMRFLib_idxval_tp *pAAt = preopt->pAAt_idxval[i]; \
						double lat = latent[i];	\
						for(int k = 0; k < pAAt->n; k++) { \
							pred[pAAt->idx[k]] += lat * pAAt->val[k]; \
						}			\
					}				\
				}					\
			}

			RUN_CODE_BLOCK(GMRFLib_MAX_THREADS(), 0, 0);
#undef CODE_BLOCK
		} else {
			// one loop

			for (int i = 0; i < preopt->n; i++) {
				GMRFLib_idxval_tp *At = preopt->At_idxval[i];
				double lat = latent[i];
				for (int k = 0; k < At->n; k++) {
					pred[offset + At->idx[k]] += lat * At->val[k];
				}
			}
		}

	} else {

		// not data-rich case

		if (preopt->pA_idxval) {
			// both loops
#define CODE_BLOCK							\
			for(int j = 0; j < 2; j++) {			\
				CODE_BLOCK_SET_THREAD_ID();		\
				if (j == 0) {				\
					for (int i = 0; i < preopt->npred; i++) { \
						if (preopt->A_idxval[i]) { \
							GMRFLib_idxval_tp *elm = preopt->A_idxval[i]; \
							DOT_PRODUCT_SERIAL(pred[offset + i], elm, latent); \
						}			\
					}				\
				} else {				\
					for (int i = 0; i < preopt->mpred; i++) { \
						if (preopt->pAA_idxval[i]) { \
							GMRFLib_idxval_tp *elm = preopt->pAA_idxval[i]; \
							DOT_PRODUCT_SERIAL(pred[i], elm, latent); \
						}			\
					}				\
				}					\
			}

			RUN_CODE_BLOCK(GMRFLib_MAX_THREADS(), 0, 0);
#undef CODE_BLOCK

		} else {
			// one loop
#define CODE_BLOCK							\
			for (int i = 0; i < preopt->npred; i++) {	\
				CODE_BLOCK_SET_THREAD_ID();		\
				if (preopt->A_idxval[i]) {		\
					GMRFLib_idxval_tp *elm = preopt->A_idxval[i]; \
					DOT_PRODUCT_SERIAL(pred[offset + i], elm, latent); \
				}					\
			}

			RUN_CODE_BLOCK(GMRFLib_MAX_THREADS(), 0, 0);
#undef CODE_BLOCK
		}
	}

	if (likelihood_only) {
		Memcpy(predictor, pred, preopt->Npred * sizeof(double));
	} else {
		Memcpy(predictor, pred, preopt->mnpred * sizeof(double));
	}
	Free(pred);

	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}

int GMRFLib_preopt_predictor_core_ORIG(double *predictor, double *latent, GMRFLib_preopt_tp * preopt, int likelihood_only)
{
	// if likelihood_only, only compute the part that is needed for the likelihood.

	// if !likelihood_only, compute the whole predictor

	GMRFLib_ENTER_ROUTINE;

	int offset = 0;
	double *pred = Calloc(preopt->mnpred, double);

	if (preopt->pA_idxval) {
		offset = preopt->mpred;
		assert(preopt->mpred > 0);
	} else {
		assert(preopt->mpred == 0);
	}

#define CODE_BLOCK							\
	for (int i = 0; i < preopt->npred; i++) {			\
		CODE_BLOCK_SET_THREAD_ID();				\
		if (preopt->A_idxval[i]) {				\
			GMRFLib_idxval_tp *elm = preopt->A_idxval[i];	\
			for (int jj = 0; jj < preopt->A_idxval[i]->n; jj++) { \
				pred[offset + i] += elm->val[jj] * latent[elm->idx[jj]]; \
			}						\
		}							\
	}

	RUN_CODE_BLOCK(GMRFLib_MAX_THREADS(), 0, 0);
#undef CODE_BLOCK

	if (preopt->mpred) {
#define CODE_BLOCK							\
		for (int i = 0; i < preopt->mpred; i++) {		\
			CODE_BLOCK_SET_THREAD_ID();			\
			if (preopt->pA_idxval[i]) {			\
				GMRFLib_idxval_tp *elm = preopt->pA_idxval[i]; \
				for (int jj = 0; jj < preopt->pA_idxval[i]->n; jj++) { \
					pred[i] += elm->val[jj] * pred[offset + elm->idx[jj]]; \
				}					\
			}						\
		}

		RUN_CODE_BLOCK(GMRFLib_MAX_THREADS(), 0, 0);
#undef CODE_BLOCK
	}

	if (likelihood_only) {
		Memcpy(predictor, pred, preopt->Npred * sizeof(double));
	} else {
		Memcpy(predictor, pred, preopt->mnpred * sizeof(double));
	}
	Free(pred);
	GMRFLib_LEAVE_ROUTINE;

	return GMRFLib_SUCCESS;
}

int GMRFLib_preopt_predictor_moments(double *mean, double *variance, GMRFLib_preopt_tp * preopt,
				     GMRFLib_problem_tp * problem, double *optional_mean)
{
	GMRFLib_ENTER_ROUTINE;

	// compute the marginal mean and variance for the linear predictor
	// either 'mean' and/or 'variance' could be NULL
	int npred = preopt->npred;
	int mpred = preopt->mpred;
	int mnpred = preopt->mnpred;
	int offset = mpred;
	int compute_mean = (mean ? 1 : 0);
	int compute_variance = (variance ? 1 : 0);
	int data_rich_case = (IMAX(preopt->mpred, preopt->npred) > preopt->n);
	double *mm = (optional_mean ? optional_mean : problem->sub_mean_constr);

	if (compute_mean) {
		Memset((void *) mean, 0, (size_t) mnpred * sizeof(double));
	}
	if (compute_variance) {
		Memset((void *) variance, 0, (size_t) mnpred * sizeof(double));
	}

	if (!compute_mean && !compute_variance) {
		GMRFLib_LEAVE_ROUTINE;
		return GMRFLib_SUCCESS;
	}

	if (preopt->pA) {
		if (compute_mean && !compute_variance) {

			// mean only
			if (data_rich_case) {
				for (int i = 0; i < preopt->n; i++) {
					GMRFLib_idxval_tp *pAAt = preopt->pAAt_idxval[i];
					double lat = mm[i];
					for (int k = 0; k < pAAt->n; k++) {
						mean[pAAt->idx[k]] += lat * pAAt->val[k];
					}
				}
			} else {
#define CODE_BLOCK							\
				for(int i = 0; i < mpred; i++) {	\
					CODE_BLOCK_SET_THREAD_ID();	\
					GMRFLib_idxval_tp *elm = preopt->pAA_idxval[i]; \
					DOT_PRODUCT_SERIAL(mean[i], elm, mm); \
				}

				RUN_CODE_BLOCK(GMRFLib_MAX_THREADS(), 0, 0);
#undef CODE_BLOCK
			}
		} else {
			// both mean and variance

#define CODE_BLOCK							\
			for(int i = 0; i < mpred; i++) {		\
				CODE_BLOCK_SET_THREAD_ID();		\
				double var = 0.0, *cov;			\
				int k, j, kk, jj;			\
				GMRFLib_idxval_tp *elm = preopt->pAA_idxval[i]; \
				for(k = 0; k < preopt->pAA_idxval[i]->n; k++) {	\
					j = elm->idx[k];		\
					if (compute_mean) {		\
						mean[i] += elm->val[k] * mm[j];	\
					}				\
					cov = GMRFLib_Qinv_get(problem, j, j); \
					var += SQR(elm->val[k]) * *cov;	\
					double tvar = 0.0;		\
					for(kk = k+1; kk < preopt->pAA_idxval[i]->n; kk++){ \
						jj = elm->idx[kk];	\
						cov = GMRFLib_Qinv_get(problem, j, jj);	\
						tvar += elm->val[kk] * *cov; \
					}				\
					var += 2.0 * elm->val[k] * tvar; \
				}					\
				variance[i] = var;			\
			}

			RUN_CODE_BLOCK(GMRFLib_MAX_THREADS(), 0, 0);
#undef CODE_BLOCK
		}
	}

	int err_count = 0;
	if (compute_mean && !compute_variance) {

		// mean only
		if (data_rich_case) {
			for (int i = 0; i < preopt->n; i++) {
				GMRFLib_idxval_tp *At = preopt->At_idxval[i];
				double lat = mm[i];
				for (int k = 0; k < At->n; k++) {
					mean[offset + At->idx[k]] += lat * At->val[k];
				}
			}
		} else {
#define CODE_BLOCK							\
			for(int i = 0; i < npred; i++) {		\
				CODE_BLOCK_SET_THREAD_ID();		\
				GMRFLib_idxval_tp *elm = preopt->A_idxval[i]; \
				for(int k = 0; k < preopt->A_idxval[i]->n; k++){ \
					mean[offset + i] += elm->val[k] * mm[elm->idx[k]]; \
				}					\
			}

			RUN_CODE_BLOCK(GMRFLib_MAX_THREADS(), 0, 0);
#undef CODE_BLOCK
		}

	} else {

		// both mean and variance
#define CODE_BLOCK							\
		for(int i = 0; i < npred; i++) {			\
			CODE_BLOCK_SET_THREAD_ID();			\
			double var = 0.0, zero = 0.0, *cov = NULL;	\
			int k, j, kk, jj;				\
			GMRFLib_idxval_tp *elm = preopt->A_idxval[i];	\
			for(k = 0; k < preopt->A_idxval[i]->n; k++){	\
				j = elm->idx[k];			\
				if (compute_mean) {			\
					mean[offset + i] += elm->val[k] * mm[j]; \
				}					\
				cov = GMRFLib_Qinv_get(problem, j, j);	\
				var += SQR(elm->val[k]) * *cov;		\
				double tvar = 0.0;			\
				for(kk = k+1; kk < preopt->A_idxval[i]->n; kk++){ \
					jj = elm->idx[kk];		\
					cov = GMRFLib_Qinv_get(problem, j, jj);	\
					if (!cov) {			\
						err_count++;		\
						cov = &zero;		\
					}				\
					tvar += elm->val[kk] * *cov;	\
				}					\
				var += 2.0 * elm->val[k] * tvar;	\
			}						\
			variance[offset + i] = var;			\
		}

		RUN_CODE_BLOCK(GMRFLib_MAX_THREADS(), 0, 0);
#undef CODE_BLOCK
	}

	if (err_count) {
		char *msg = NULL;
		GMRFLib_sprintf(&msg, "Missing %1d covariances; The A-matrix has not the proper rank. Please check.", err_count);
		GMRFLib_ERROR_MSG(GMRFLib_EMISC, msg);
	}

	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}

int GMRFLib_preopt_predictor_moments_ORIG(double *mean, double *variance, GMRFLib_preopt_tp * preopt,
					  GMRFLib_problem_tp * problem, double *optional_mean)
{
	GMRFLib_ENTER_ROUTINE;

	// compute the marginal mean and variance for the linear predictor
	// either 'mean' and/or 'variance' could be NULL
	int npred = preopt->npred;
	int mpred = preopt->mpred;
	int mnpred = preopt->mnpred;
	int offset = mpred;
	int compute_mean = (mean ? 1 : 0);
	int compute_variance = (variance ? 1 : 0);

	double *mm = (optional_mean ? optional_mean : problem->sub_mean_constr);

	if (compute_mean) {
		Memset((void *) mean, 0, (size_t) mnpred * sizeof(double));
	}
	if (compute_variance) {
		Memset((void *) variance, 0, (size_t) mnpred * sizeof(double));
	}

	if (!compute_mean && !compute_variance) {
		GMRFLib_LEAVE_ROUTINE;
		return GMRFLib_SUCCESS;
	}

	if (compute_mean && !compute_variance) {
#define CODE_BLOCK							\
		for(int i = 0; i < mpred; i++) {			\
			CODE_BLOCK_SET_THREAD_ID();			\
			GMRFLib_idxval_tp *elm = preopt->pAA_idxval[i]; \
			for(int k = 0; k < preopt->pAA_idxval[i]->n; k++) { \
				mean[i] += elm->val[k] * mm[elm->idx[k]]; \
			}						\
		}

		RUN_CODE_BLOCK(GMRFLib_MAX_THREADS(), 0, 0);
#undef CODE_BLOCK

	} else {

#define CODE_BLOCK							\
		for(int i = 0; i < mpred; i++) {			\
			CODE_BLOCK_SET_THREAD_ID();			\
			double var = 0.0, *cov;				\
			int k, j, kk, jj;				\
			GMRFLib_idxval_tp *elm = preopt->pAA_idxval[i]; \
			for(k = 0; k < preopt->pAA_idxval[i]->n; k++) {	\
				j = elm->idx[k];			\
				if (compute_mean) {			\
					mean[i] += elm->val[k] * mm[j];	\
				}					\
				cov = GMRFLib_Qinv_get(problem, j, j);	\
				var += SQR(elm->val[k]) * *cov;		\
				double tvar = 0.0;			\
				for(kk = k+1; kk < preopt->pAA_idxval[i]->n; kk++){ \
					jj = elm->idx[kk];		\
					cov = GMRFLib_Qinv_get(problem, j, jj);	\
					tvar += elm->val[kk] * *cov;	\
				}					\
				var += 2.0 * elm->val[k] * tvar;	\
			}						\
			variance[i] = var;				\
		}

		RUN_CODE_BLOCK(GMRFLib_MAX_THREADS(), 0, 0);
#undef CODE_BLOCK
	}

	int err_count = 0;
	if (compute_mean && !compute_variance) {
#define CODE_BLOCK							\
		for(int i = 0; i < npred; i++) {			\
			CODE_BLOCK_SET_THREAD_ID();			\
			GMRFLib_idxval_tp *elm = preopt->A_idxval[i];	\
			for(int k = 0; k < preopt->A_idxval[i]->n; k++){ \
				mean[offset + i] += elm->val[k] * mm[elm->idx[k]]; \
			}						\
		}

		RUN_CODE_BLOCK(GMRFLib_MAX_THREADS(), 0, 0);
#undef CODE_BLOCK
	} else {

#define CODE_BLOCK							\
		for(int i = 0; i < npred; i++) {			\
			CODE_BLOCK_SET_THREAD_ID();			\
			double var = 0.0, zero = 0.0, *cov = NULL;	\
			int k, j, kk, jj;				\
			GMRFLib_idxval_tp *elm = preopt->A_idxval[i];	\
			for(k = 0; k < preopt->A_idxval[i]->n; k++){	\
				j = elm->idx[k];				\
				if (compute_mean) {			\
					mean[offset + i] += elm->val[k] * mm[j]; \
				}					\
				cov = GMRFLib_Qinv_get(problem, j, j);	\
				var += SQR(elm->val[k]) * *cov;		\
				double tvar = 0.0;			\
				for(kk = k+1; kk < preopt->A_idxval[i]->n; kk++){ \
					jj = elm->idx[kk];		\
					cov = GMRFLib_Qinv_get(problem, j, jj);	\
					if (!cov) {			\
						err_count++;		\
						cov = &zero;		\
					}				\
					tvar += elm->val[kk] * *cov;	\
				}					\
				var += 2.0 * elm->val[k] * tvar;	\
			}						\
			variance[offset + i] = var;			\
		}

		RUN_CODE_BLOCK(GMRFLib_MAX_THREADS(), 0, 0);
#undef CODE_BLOCK
	}

	if (err_count) {
		char *msg = NULL;
		GMRFLib_sprintf(&msg, "Missing %1d covariances; The A-matrix has not the proper rank. Please check.", err_count);
		GMRFLib_ERROR_MSG(GMRFLib_EMISC, msg);
	}

	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}

int GMRFLib_preopt_update(GMRFLib_preopt_tp * preopt, double *like_b, double *like_c)
{
	int id = GMRFLib_thread_id, np = preopt->Npred;

	if (!(preopt->like_b[id])) {
		preopt->like_b[id] = Calloc(np, double);
		preopt->like_c[id] = Calloc(np, double);
	}
	Memcpy(preopt->like_b[id], like_b, np * sizeof(double));
	Memcpy(preopt->like_c[id], like_c, np * sizeof(double));

	if (!(preopt->total_b[id])) {
		preopt->total_b[id] = Calloc(preopt->n, double);
	}
	Memset(preopt->total_b[id], 0, preopt->n * sizeof(double));
	GMRFLib_preopt_bnew(preopt->total_b[id], preopt);

	return GMRFLib_SUCCESS;
}

int GMRFLib_preopt_free(GMRFLib_preopt_tp * preopt)
{
	int id = GMRFLib_thread_id;

	if (!preopt) {
		return GMRFLib_SUCCESS;
	}
#pragma omp parallel sections
	{
#pragma omp section
		{
			GMRFLib_thread_id = id;
			if (preopt->pAA_idxval) {
				for (int i = 0; i < preopt->mpred; i++) {
					GMRFLib_idxval_free(preopt->pAA_idxval[i]);
				}
				Free(preopt->pAA_idxval);
			}
			for (int i = 0; i < preopt->n; i++) {
				GMRFLib_idxval_free(preopt->AtA_idxval[i][0]);
				for (int jj = 0; jj < preopt->like_graph->lnnbs[i]; jj++) {
					GMRFLib_idxval_free(preopt->AtA_idxval[i][1 + jj]);
				}
			}
			Free(preopt->AtA_idxval);
			if (preopt->pA_idxval) {
				for (int i = 0; i < preopt->mpred; i++) {
					GMRFLib_idxval_free(preopt->pA_idxval[i]);
				}
				for (int i = 0; i < preopt->n; i++) {
					GMRFLib_idxval_free(preopt->pAAt_idxval[i]);
				}
			}
			for (int i = 0; i < preopt->npred; i++) {
				GMRFLib_idxval_free(preopt->A_idxval[i]);
			}
			for (int i = 0; i < preopt->n; i++) {
				GMRFLib_idxval_free(preopt->At_idxval[i]);
			}
		}

#pragma omp section
		{
			GMRFLib_thread_id = id;
			GMRFLib_matrix_free(preopt->A);
			GMRFLib_matrix_free(preopt->pA);

			Free(preopt->idx_map_f);
			Free(preopt->idx_map_beta);
			Free(preopt->what_type);

			for (int i = 0; i < GMRFLib_MAX_THREADS(); i++) {
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
		}
	}
	Free(preopt);

	return GMRFLib_SUCCESS;
}

double *GMRFLib_preopt_measure_time(GMRFLib_preopt_tp * preopt)
{
	// return alloc'ed double *cpu measurements.
	// cpu[0] and cpu[1] is the time for doing Q %*% x.
	// cpu[0] does elementwise computations, cpu[1] does blockwise using tabulated Q.
	
	double value = 0.0;
	double *cpu = Calloc(2, double);

	GMRFLib_Qfunc_tp *like_Qfunc = preopt->like_Qfunc;
	GMRFLib_graph_tp *like_graph = preopt->like_graph;
	void *like_Qfunc_arg = preopt->like_Qfunc_arg;

	// this will be measure with serial or with group
	cpu[0] = -GMRFLib_cpu();
	for (int i = 0; i < like_graph->n; i++) {
		value += like_Qfunc(i, i, NULL, like_Qfunc_arg);
		for (int jj = 0, j; jj < like_graph->lnnbs[i]; jj++) {
			j = like_graph->lnbs[i][jj];
			value += like_Qfunc(i, j, NULL, like_Qfunc_arg);
		}
	}
	cpu[0] += GMRFLib_cpu();
	assert(!ISNAN(value));

	GMRFLib_Qfunc_tp *Qfunc = preopt->preopt_Qfunc;
	GMRFLib_graph_tp *graph = preopt->preopt_graph;
	void *Qfunc_arg = preopt->preopt_Qfunc_arg;

	Calloc_init(2 * graph->n);
	assert(calloc_work_);
	double *x = Calloc_get(graph->n);
	double *xx = Calloc_get(graph->n);
	for (int i = 0; i < graph->n; i++) {
		x[i] = GMRFLib_uniform();
	}

	GMRFLib_tabulate_Qfunc_tp *tab = NULL;
	GMRFLib_tabulate_Qfunc(&tab, graph, Qfunc, Qfunc_arg, NULL);

	// this will be measured with serial or parallel
	cpu[1] = -GMRFLib_cpu();
	GMRFLib_Qx(xx, x, graph, tab->Qfunc, tab->Qfunc_arg);
	cpu[1] += GMRFLib_cpu();

	Calloc_free();
	GMRFLib_free_tabulate_Qfunc(tab);

	return cpu;
}

double *GMRFLib_preopt_measure_time2(GMRFLib_preopt_tp * preopt)
{
	// cpu[0] measure predictor calculations for !data_rich and data_rich case
	double value = 0.0;
	double *cpu = Calloc(1, double);

	double *pred = Calloc(2 * preopt->mnpred, double);
	double *lat = pred + preopt->mnpred;

	for(int i = 0; i < preopt->mnpred; i++) {
		lat[i] = GMRFLib_uniform();
	}
	cpu[0] = -GMRFLib_cpu();
	GMRFLib_preopt_full_predictor(pred, lat, preopt);
	cpu[0] += GMRFLib_cpu();
	Free(pred);

	return cpu;
}











	// register those matrices that are empty
	static map_vpi cache;
	static int cache_must_init = 1;

	if (cache_must_init) {
#pragma omp critical (Name_0da3ba144656768918b2515c57b703c3d8e57f2a)
		if (cache_must_init) {
			map_vpi_init_hint(&cache, 128);
			cache_must_init = 0;
		}
	}

#define CODE_BLOCK						\
	for (int i = 0; i < preopt->n; i++) {			\
		if (A[i]) {					\
			GMRFLib_idxval_tp *elm = A[i];		\
			DOT_PRODUCT(bnew[i], elm, blike);	\
		}						\
	}

	int *ival = map_vpi_ptr(&cache, (void *) A);

	if (ival) P(*ival);
	if (!ival) FIXME("NEED TO BE CHECKED");
	
	if (!ival) {
		// this matrix needs to be checked
		int c = 0;
		for (int i = 0; i < preopt->n && !c; i++) {
			if (A[i]) c = 1;
		}
		int ival_code = (c == 1 ? 2 : 1);
		P(ival_code);

#pragma omp critical (Name_2e46eda7cf12cd72251f668c56b12756c4f409d5)
		{
			map_vpi_set(&cache, (void *) A, ival_code);
		}

		if (ival_code == 1) {
			FIXME("bnew = 0");
			Memset(bnew, 0, preopt->n * sizeof(double));
		} else if (ival_code == 2) {
			RUN_CODE_BLOCK(GMRFLib_MAX_THREADS(), 0, 0);
		} else {
			assert(0 == 1);
		}
	} else if (*ival == 1) {
			FIXME("bnew = 0");
			Memset(bnew, 0, preopt->n * sizeof(double));
	} else if (*ival == 2) {
			RUN_CODE_BLOCK(GMRFLib_MAX_THREADS(), 0, 0);
	} else {
		assert(0 == 1);
	}
		

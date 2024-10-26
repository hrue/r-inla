
/* pre-opt.c
 * 
 * Copyright (C) 2021-2024 Havard Rue
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more detail.
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

#include <time.h>
#include <strings.h>
#include <stdlib.h>
#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

#if !defined(WINDOWS)
#include <unistd.h>
#include <sys/types.h>
#include <pwd.h>
#endif

int GMRFLib_preopt_init(GMRFLib_preopt_tp **preopt, int npred, int nf, int **c, double **w,
			GMRFLib_graph_tp **f_graph, GMRFLib_Qfunc_tp **f_Qfunc, void **f_Qfunc_arg,
			char *f_sumzero, GMRFLib_constr_tp **f_constr, double *f_diag,
			GMRFLib_Qfunc_tp ***ff_Qfunc, void ***ff_Qfunc_arg, GMRFLib_matrix_tp **f_Alocal,
			int nbeta, double **covariate,
			double *prior_precision, GMRFLib_bfunc_tp **bfunc, GMRFLib_ai_param_tp *UNUSED(ai_par),
			char *pA_fnm, GMRFLib_matrix_tp **global_constr)
{
	assert(omp_get_thread_num() == 0);

#define SHOW_TIME(_msg)							\
	if (debug) {							\
		printf("\t\tGMRFLib_preopt_init: %-16s %7.2fs\n", _msg, GMRFLib_timer() - tref); \
		tref =  GMRFLib_timer();					\
	}

	if (!preopt) {
		return GMRFLib_SUCCESS;
	}
	GMRFLib_ENTER_ROUTINE;

	int N = 0, *idx_map_f = NULL, *idx_map_beta = NULL, offset, nrow = 0, ncol = 0;
	int debug = GMRFLib_DEBUG_IF_TRUE();
	const int debug_detailed = 0;
	const int do_prune = 1;

	double tref = GMRFLib_timer();
	double **ww = NULL;
	GMRFLib_constr_tp *fc = NULL;
	GMRFLib_idx_tp **pAA_pattern = NULL;
	GMRFLib_idxval_tp **A_idxval = NULL;
	GMRFLib_idxval_tp ***AtA_idxval = NULL;
	GMRFLib_idxval_tp **At_idxval = NULL;
	GMRFLib_idxval_tp **pAA_idxval = NULL;
	GMRFLib_idxval_tp **pAAt_idxval = NULL;
	GMRFLib_idxval_tp **pA_idxval = NULL;
	GMRFLib_idxval_tp *elm = NULL;
	GMRFLib_matrix_tp *pA = NULL;

	ww = Calloc(nf, double *);
	for (int i = 0; i < nf; i++) {
		ww[i] = Calloc(npred, double);
		GMRFLib_fill(npred, 1.0, ww[i]);
	}
	if (w) {
		for (int i = 0; i < nf; i++) {
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
		for (int i = 0; i < nf; i++) {
			idx_map_f[i] = offset;
			offset += f_graph[i]->n;
		}
		idx_map_f[nf] = offset;
	}
	if (nbeta) {
		idx_map_beta = Calloc(nbeta + 1, int);
		for (int i = 0; i < nbeta; i++) {
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
		for (int j = 0; j < nf; j++) {
			GMRFLib_ged_insert_graph(ged, f_graph[j], idx_map_f[j]);
		}
	}
	for (int i = 0; i < nf; i++) {
		for (int j = i + 1; j < nf; j++) {
			if (ff_Qfunc && ff_Qfunc[i][j]) {
				for (int k = 0; k < IMIN(f_graph[i]->n, f_graph[j]->n); k++) {
					GMRFLib_ged_add(ged, idx_map_f[i] + k, idx_map_f[j] + k);
				}
			}
		}
	}
	if (nbeta) {
		for (int i = 0; i < nbeta; i++) {
			GMRFLib_ged_add(ged, idx_map_beta[i], idx_map_beta[i]);
		}
	}

	GMRFLib_ged_build(&((*preopt)->latent_graph), ged);
	GMRFLib_ged_free(ged);
	// not needed as its only one option
	// (*preopt)->latent_Qfunc = GMRFLib_preopt_latent_Qfunc;
	(*preopt)->latent_Qfunc_arg = (void *) *preopt;

	/*
	 * build the constraint, if any. 
	 */
	int nconstr = 0;
	if (nf && f_sumzero) {
		for (int k = 0; k < nf; k++) {
			nconstr += (f_sumzero[k] ? 1 : 0);
		}
	}

	if (nf && f_constr) {
		for (int k = 0; k < nf; k++) {
			fc = f_constr[k];
			if (fc) {
				nconstr += fc->nc;
			}
		}
	}

	int ngc = 0;
	if (global_constr) {
		ngc = global_constr[1]->nrow;
		nconstr += ngc;
	}

	if (nconstr) {
		int constr_no = 0;
		int nn = (*preopt)->latent_graph->n;
		GMRFLib_constr_tp *constr = NULL;

		GMRFLib_make_empty_constr(&constr);
		constr->a_matrix = Calloc(nn * nconstr, double);
		constr->e_vector = Calloc(nconstr, double);
		constr->nc = nconstr;

		if (nf && f_sumzero) {
			for (int k = 0; k < nf; k++) {
				if (f_sumzero[k]) {
					for (int i = idx_map_f[k]; i < idx_map_f[k + 1]; i++) {
						constr->a_matrix[i * constr->nc + constr_no] = 1.0;
					}
					constr->e_vector[constr_no] = 0.0;
					constr_no++;
				}
			}
		}
		if (nf && f_constr) {
			for (int k = 0; k < nf; k++) {
				fc = f_constr[k];
				if (fc) {
					for (int j = 0; j < fc->nc; j++) {
						for (int i = idx_map_f[k], ii = 0; i < idx_map_f[k + 1]; i++, ii++) {
							constr->a_matrix[i * constr->nc + constr_no] = fc->a_matrix[ii * fc->nc + j];
						}
						constr->e_vector[constr_no] = fc->e_vector[j];
						constr_no++;
					}
				}
			}
		}

		if (ngc) {
			// validate. we cannot do this before now
			if (!(ngc * nn == global_constr[0]->nrow)) {
				fprintf(stderr, "\n\n\n");
				fprintf(stderr, "\t*** Number of global constraints = %1d,  but the length = %1d\n", ngc, global_constr[0]->nrow);
				fprintf(stderr, "\t*** is not a multiplum of the size of the latent = %1d, exit...\n", nn);
				GMRFLib_ASSERT(ngc * nn == global_constr[0]->nrow, GMRFLib_ESNH);
			}

			for (int j = 0; j < ngc; j++) {
				int off = j * nn;
				for (int i = 0; i < nn; i++) {
					constr->a_matrix[i * constr->nc + constr_no] = global_constr[0]->A[i + off];
				}
				constr->e_vector[constr_no] = global_constr[1]->A[j];
				constr_no++;
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
	for (int i = 0; i < N; i++) {
		(*preopt)->what_type[i] = GMRFLib_preopt_what_type(i, *preopt);
	}

	if (debug_detailed) {
		printf("\tnpred %1d nf %1d nbeta %1d\n", npred, nf, nbeta);
		for (int i = 0; i < npred; i++) {
			printf("data %1d\n", i);
			for (int j = 0; j < nf; j++) {
				printf("\t\tf[%1d]  index %1d  weight %.6f\n", j, c[j][i], ww[j][i]);
			}
			for (int j = 0; j < nbeta; j++) {
				printf("\t\tbeta[%1d]  x %.6f\n", j, covariate[j][i]);
			}
		}
	}

	SHOW_TIME("setup2");

	// build up structure for the likelihood part

	GMRFLib_ASSERT(GMRFLib_OPENMP_IN_PARALLEL_ONEPLUS_THREAD() == 0, GMRFLib_ESNH);
	A_idxval = GMRFLib_idxval_ncreate(npred);

#pragma omp parallel for num_threads(GMRFLib_openmp->max_threads_outer)
	for (int i = 0; i < npred; i++) {
		GMRFLib_idxval_tp *iv = NULL;
		double val = 0.0;
		int idx = 0;

		for (int jj = 0; jj < nf; jj++) {
			int have_Alocal = ((f_Alocal != NULL) && (f_Alocal[jj] != NULL) && LEGAL(i, f_Alocal[jj]->nrow));
			int nn = f_graph[jj]->n;
			int off = idx_map_f[jj];

			val = ww[jj][i];
			if (c[jj][i] >= 0 && !ISZERO(val)) {
				idx = c[jj][i] + off;
				GMRFLib_idxval_add(&(A_idxval[i]), idx, val);
			}

			if (have_Alocal) {
				if (i == 0) {
					// turn this one off. not sure its a good idea
					if (0 && (nn != f_Alocal[jj]->nrow)) {
						fprintf(stderr, "\n\n%s:%1d (%s) warning: (nrow(Alocal) = %1d) != (f[%1d]->n = %1d)\n",
							__FILE__, __LINE__, __GMRFLib_FuncName, f_Alocal[jj]->nrow, jj, nn);
						fprintf(stderr, "\t\tThis might, or might not, be what you want\n\n");
						fflush(stderr);
					}
				}

				if (iv) {
					iv->n = 0;
				}
				GMRFLib_matrix_get_row_idxval(&iv, i, f_Alocal[jj], 0);
				if (iv) {
					if (iv->n && !LEGAL(iv->idx[iv->n - 1], nn)) {
						fprintf(stderr, "%s:%1d (%s) error: Illegal index in Alocal: idx.f = %1d idx = %1d n = %1d\n",
							__FILE__, __LINE__, __GMRFLib_FuncName, jj, iv->idx[iv->n - 1], nn);
						assert(LEGAL(iv->idx[iv->n - 1], nn));
					}

					for (int ii = 0; ii < iv->n; ii++) {
						val = iv->val[ii];
						if (!ISZERO(val)) {
							idx = iv->idx[ii] + off;
							GMRFLib_idxval_add(&(A_idxval[i]), idx, val);
						}
					}
				}
			}
		}
		GMRFLib_idxval_free(iv);

		for (int jj = 0; jj < nbeta; jj++) {
			val = covariate[jj][i];
			if (!ISZERO(val)) {
				idx = idx_map_beta[jj];
				GMRFLib_idxval_add(&(A_idxval[i]), idx, val);
			}
		}
		GMRFLib_idxval_prepare(&(A_idxval[i]), 1, 1);
		if (do_prune) {
			GMRFLib_idxval_prune(A_idxval[i]);
		}
	}

	GMRFLib_idxval_to_matrix(&((*preopt)->A), A_idxval, npred, N);
	SHOW_TIME("A_idxval");

	// need also At_.. below, if (pA)
	At_idxval = GMRFLib_idxval_ncreate(N);
	for (int i = 0; i < npred; i++) {
		elm = A_idxval[i];
		for (int k = 0; k < A_idxval[i]->n; k++) {
			GMRFLib_idxval_add(&(At_idxval[elm->idx[k]]), i, elm->val[k]);
		}
	}
	GMRFLib_idxval_prepare(At_idxval, N, GMRFLib_MAX_THREADS());
	if (do_prune) {
		GMRFLib_idxval_nprune(At_idxval, N, GMRFLib_MAX_THREADS());
	}

	SHOW_TIME("At_idxval");
	if (debug_detailed) {
		for (int i = 0; i < npred; i++) {
			GMRFLib_idxval_printf(stdout, A_idxval[i], "A_idxval");
		}
		for (int i = 0; i < N; i++) {
			GMRFLib_idxval_printf(stdout, At_idxval[i], "At_idxval");
		}
	}

	if (pA_fnm) {
		pA = GMRFLib_read_fmesher_file(pA_fnm, (long int) 0, -1);
		assert(pA);

		if (debug_detailed) {
			printf("read pA from [%s]\n", pA_fnm);
			printf("\tnrow %d ncol %d nelms %d\n", pA->nrow, pA->ncol, pA->elems);
			for (int i = 0; i < pA->elems; i++) {
				printf("\ti j x %d %d %f\n", pA->i[i], pA->j[i], pA->values[i]);
			}
		}

		nrow = pA->nrow;
		ncol = pA->ncol;
		assert(ncol == npred);
		SHOW_TIME("read pA");

		// this is need to compute the linear predictor later
		pA_idxval = GMRFLib_idxval_ncreate(nrow);
		for (int k = 0; k < pA->elems; k++) {
			int i = pA->i[k];
			int j = pA->j[k];
			GMRFLib_idxval_add(&(pA_idxval[i]), j, pA->values[k]);
		}
		GMRFLib_idxval_prepare(pA_idxval, nrow, GMRFLib_MAX_THREADS());
		if (do_prune) {
			GMRFLib_idxval_nprune(pA_idxval, nrow, GMRFLib_MAX_THREADS());
		}
		(*preopt)->pA = pA;
		SHOW_TIME("create pA_idxval");

		// to avoid to much 'realloc', I can compute the the length in 'm' and then add terms
		pAA_pattern = Calloc(nrow, GMRFLib_idx_tp *);

		// this will keep the working 'idxval' within the thread, and we can free it at the end
		GMRFLib_idxval_tp **row_idxval_hold = Calloc(GMRFLib_MAX_THREADS(), GMRFLib_idxval_tp *);

#pragma omp parallel for num_threads(GMRFLib_openmp->max_threads_outer)
		for (int i = 0; i < nrow; i++) {
			int thread = omp_get_thread_num();
			GMRFLib_idxval_tp *row_idxval = row_idxval_hold[thread];
			if (row_idxval) {
				// we do not free it, we can just pretend its empty and use it again
				row_idxval->n = 0;
			}
			// last argument = 0, as we do not need to sort it
			GMRFLib_matrix_get_row_idxval(&row_idxval, i, pA, 0);

			// total length
			int m = 0;
			for (int jj = 0; jj < row_idxval->n; jj++) {
				int j = row_idxval->idx[jj];
				m += A_idxval[j]->n;
			}
			GMRFLib_idx_create_x(&(pAA_pattern[i]), m);

			for (int jj = 0; jj < row_idxval->n; jj++) {
				int j = row_idxval->idx[jj];
				// use the _nadd to append a whole vector
				GMRFLib_idx_nadd(&(pAA_pattern[i]), A_idxval[j]->n, A_idxval[j]->idx);
				// instead of this old code
				// for (int kk = 0; kk < A_idxval[j]->n; kk++) {
				// int k = A_idxval[j]->idx[kk];
				// GMRFLib_idx_add(&(pAA_pattern[i]), k); }
			}
			GMRFLib_idx_uniq(pAA_pattern[i]);      /* this also sorts */
			GMRFLib_idxval_free(row_idxval);
		}

		for (int i = 0; i < GMRFLib_MAX_THREADS(); i++) {
			if (row_idxval_hold[i]) {
				GMRFLib_idxval_free(row_idxval_hold[i]);
			}
		}
		Free(row_idxval_hold);
		SHOW_TIME("pAA_pattern");


		if (debug_detailed) {
			char *crow = Calloc(N + 1, char);
			crow[N] = '\0';
			for (int i = 0; i < nrow; i++) {
				Memset(crow, ' ', N * sizeof(char));
				for (int k = 0; k < pAA_pattern[i]->n; k++) {
					int j = pAA_pattern[i]->idx[k];
					// printf("Add crow i j %d %d\n", i, j);
					crow[j] = '.';
				}
				printf("pAA%2d [%s]\n", i, crow);
			}
			Free(crow);
		}
		// first make a empty one filled with zeros to get the pattern. since pAA_pattern is sorted, then this will be sorted as well
		pAA_idxval = GMRFLib_idxval_ncreate(nrow);

#pragma omp parallel for num_threads(GMRFLib_openmp->max_threads_outer)
		for (int i = 0; i < nrow; i++) {
			int *idx = pAA_pattern[i]->idx;
			for (int k = 0; k < pAA_pattern[i]->n; k++) {
				GMRFLib_idxval_add(&(pAA_idxval[i]), idx[k], 0.0);
			}
		}

		SHOW_TIME("init pAA_idxval");

#pragma omp parallel for num_threads(GMRFLib_openmp->max_threads_outer)
		for (int i = 0; i < nrow; i++) {
			int step;
			int steps[] = { 262144, 32768, 4096, 512, 64, 8, 1 };
			int nsteps = sizeof(steps) / sizeof(int);
			int row_n;
			GMRFLib_idxval_tp *row_idxval = NULL;
			GMRFLib_idxval_tp *row_elm = NULL;

			GMRFLib_matrix_get_row_idxval(&row_idxval, i, pA, 1);
			row_elm = row_idxval;
			row_n = row_idxval->n;

			for (int jj = 0; jj < pAA_pattern[i]->n; jj++) {
				int j = pAA_pattern[i]->idx[jj];
				GMRFLib_idxval_tp *At_elm = At_idxval[j];
				int At_n = At_idxval[j]->n;
				int irow = 0;
				int iAt = 0;
				while (irow < row_n && iAt < At_n) {
					int k = row_elm->idx[irow];
					int kk = At_elm->idx[iAt];
					if (k < kk) {
						irow++;
						for (int s = 0; s < nsteps; s++) {
							step = steps[s];
							if (step < row_n) {
								int ia = irow + step;
								while (ia < row_n && row_elm->idx[ia] < kk)
									ia += step;
								irow = ia - step;
							}
						}
					} else if (k > kk) {
						iAt++;
						for (int s = 0; s < nsteps; s++) {
							step = steps[s];
							if (step < At_n) {
								int ia = iAt + step;
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
		GMRFLib_idxval_prepare(pAA_idxval, nrow, GMRFLib_MAX_THREADS());
		if (do_prune) {
			GMRFLib_idxval_nprune(pAA_idxval, nrow, GMRFLib_MAX_THREADS());
		}
		SHOW_TIME("pAA_idxval");

		pAAt_idxval = GMRFLib_idxval_ncreate(N);
		for (int i = 0; i < nrow; i++) {
			elm = pAA_idxval[i];
			for (int k = 0; k < pAA_idxval[i]->n; k++) {
				GMRFLib_idxval_add(&(pAAt_idxval[elm->idx[k]]), i, elm->val[k]);
			}
		}
		GMRFLib_idxval_prepare(pAAt_idxval, N, GMRFLib_MAX_THREADS());
		if (do_prune) {
			GMRFLib_idxval_nprune(pAAt_idxval, N, GMRFLib_MAX_THREADS());
		}
		SHOW_TIME("pAAt_idxval");

		if (debug_detailed) {
			for (int i = 0; i < nrow; i++) {
				P(i);
				GMRFLib_idxval_printf(stdout, pAA_idxval[i], "pAA");
			}
			for (int i = 0; i < N; i++) {
				P(i);
				GMRFLib_idxval_printf(stdout, pAAt_idxval[i], "pAAt");
			}
		}

		for (int i = 0; i < nrow; i++) {
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
#pragma omp parallel for num_threads(GMRFLib_openmp->max_threads_outer)
	for (int i = 0; i < gen_len_At; i++) {
		for (int kk = 0; kk < gen_At[i]->n; kk++) {
			int k = gen_At[i]->idx[kk];
			for (int jj = 0; jj < gen_A[k]->n; jj++) {
				int j = gen_A[k]->idx[jj];
				if (j > i) {
					GMRFLib_ged_add(ged, i, j);
				}
			}
		}
	}
	SHOW_TIME("build graph part 1");

	GMRFLib_ged_build(&g, ged);
	GMRFLib_ged_free(ged);
	assert(g->n == gen_len_At);
	SHOW_TIME("build graph part 2");

	AtA_idxval = Calloc(gen_len_At, GMRFLib_idxval_tp **);
	for (int i = 0; i < g->n; i++) {
		AtA_idxval[i] = GMRFLib_idxval_ncreate(1 + g->lnnbs[i]);
	}

#pragma omp parallel for num_threads(GMRFLib_openmp->max_threads_outer)
	for (int i = 0; i < gen_len_At; i++) {
		int guess[2];
		int m = g->lnnbs[i];
		int *arr = g->lnbs[i];
		guess[0] = guess[1] = 0;

		for (int kk = 0; kk < gen_At[i]->n; kk++) {
			int k = gen_At[i]->idx[kk];
			for (int jj = 0; jj < gen_A[k]->n; jj++) {
				int j = gen_A[k]->idx[jj];
				if (j >= i) {
					int index = 0;
					if (i != j) {
						index = 1 + GMRFLib_iwhich_sorted_g2(j, arr, m, guess);
						assert(index > 0);
					}
					double value = gen_At[i]->val[kk] * gen_A[k]->val[jj];
					GMRFLib_idxval_add(&(AtA_idxval[i][index]), k, value);
				}
			}
		}
	}
	SHOW_TIME("AtA_idxval");

	if (debug_detailed) {
		FIXME("AtA");
		for (int i = 0; i < N; i++) {
			double sum = 0.0;
			printf("term %d %d\n", i, i);
			for (int kk = 0; kk < AtA_idxval[i][0]->n; kk++) {
				printf("\tkk idx val %d %d %f\n", kk, AtA_idxval[i][0]->idx[kk], AtA_idxval[i][0]->val[kk]);
				sum += AtA_idxval[i][0]->val[kk];
			}
			printf("\tsum %g\n", sum);
			for (int jj = 0; jj < g->lnnbs[i]; jj++) {
				int j = g->lnbs[i][jj];
				printf("term %d %d\n", i, j);
				sum = 0.0;
				for (int kk = 0; kk < AtA_idxval[i][1 + jj]->n; kk++) {
					printf("\tkk idx val %d %d %f\n", kk, AtA_idxval[i][1 + jj]->idx[kk], AtA_idxval[i][1 + jj]->val[kk]);
					sum += AtA_idxval[i][1 + jj]->val[kk];
				}
				printf("\tsum %g\n", sum);
			}
		}
	}
#pragma omp parallel for num_threads(GMRFLib_openmp->max_threads_outer)
	for (int i = 0; i < g->n; i++) {
		GMRFLib_idxval_prepare(AtA_idxval[i], 1 + g->lnnbs[i], 1);
		if (do_prune) {
			GMRFLib_idxval_nprune(AtA_idxval[i], 1 + g->lnnbs[i], 1);
		}
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
	// not needed as its only one option
	// (*preopt)->like_Qfunc = GMRFLib_preopt_like_Qfunc;
	// (*preopt)->like_Qfunc_k = GMRFLib_preopt_like_Qfunc_k;
	(*preopt)->bfunc = bfunc;
	(*preopt)->nf = (*preopt)->n - nbeta;
	(*preopt)->nbeta = nbeta;

	GMRFLib_graph_tp *g_arr[2];
	g_arr[0] = (*preopt)->latent_graph;
	g_arr[1] = (*preopt)->like_graph;
	GMRFLib_graph_union(&((*preopt)->preopt_graph), g_arr, 2);

#if !defined(WINDOWS)
	if (getenv("INLA_INTERNAL_DUMP_GRAPH")) {
		static int count = 0;
		char *homedir = getenv("HOME");
		if (!homedir) {
			homedir = getpwuid(getuid())->pw_dir;
		}
		if (!homedir) {
			homedir = strdup("./");
		}
		char *fnm = NULL;
		GMRFLib_sprintf(&fnm, "%s/INLA-graph-pid%1d-count%1d.txt", homedir, (int) getpid(), ++count);
		GMRFLib_graph_write(fnm, (*preopt)->preopt_graph);
		fprintf(stderr, "\n\t*** write graph to file [%s]\n", fnm);
	}
#endif

	(*preopt)->preopt_Qfunc = GMRFLib_preopt_Qfunc;
	(*preopt)->preopt_Qfunc_arg = (void *) *preopt;
	(*preopt)->gcpo_Qfunc = GMRFLib_preopt_gcpo_Qfunc;

	// tabulate _is_nb() as this will speedup _preopt_Qfunc()
	(*preopt)->preopt_graph_latent_is_nb = Calloc((*preopt)->preopt_graph->n, char *);
	(*preopt)->preopt_graph_like_is_nb = Calloc((*preopt)->preopt_graph->n, char *);

#pragma omp parallel for num_threads(GMRFLib_openmp->max_threads_outer)
	for (int i = 0; i < (*preopt)->preopt_graph->n; i++) {
		int lnnbs = (*preopt)->preopt_graph->lnnbs[i];
		if (lnnbs) {
			char *store = Calloc(2 * lnnbs, char);
			(*preopt)->preopt_graph_latent_is_nb[i] = store;
			(*preopt)->preopt_graph_like_is_nb[i] = store + lnnbs;
			for (int k = 0; k < lnnbs; k++) {
				int j = (*preopt)->preopt_graph->lnbs[i][k];
				if (GMRFLib_graph_is_nb(i, j, (*preopt)->latent_graph)) {
					(*preopt)->preopt_graph_latent_is_nb[i][k] = 1;
				}
				if (GMRFLib_graph_is_nb(i, j, (*preopt)->like_graph)) {
					(*preopt)->preopt_graph_like_is_nb[i][k] = 1;
				}
			}
		}
	}

	for (int i = 0; i < nf; i++) {
		Free(ww[i]);
	}
	Free(ww);

	GMRFLib_openmp_implement_strategy(GMRFLib_OPENMP_PLACES_DEFAULT, NULL, NULL);

	SHOW_TIME("admin2");
#undef  SHOW_TIME

	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}

GMRFLib_preopt_type_tp GMRFLib_preopt_what_type(int node, GMRFLib_preopt_tp *preopt)
{
	GMRFLib_preopt_type_tp t = { GMRFLib_PREOPT_TP___VOID, -1, -1 };
	if (preopt->nf && node < preopt->idx_map_f[preopt->nf]) {
		t.tp = GMRFLib_PREOPT_TP_F;
		for (int i = 0; i < preopt->nf; i++) {
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

double GMRFLib_preopt_latent_Qfunc(int thread_id, int node, int nnode, double *UNUSED(values), void *arg)
{
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
	{
		switch (jt.tp) {
		case GMRFLib_PREOPT_TP_F:
		{
			if (same_tp) {
				if (same_idx || GMRFLib_graph_is_nb(it.idx, jt.idx, a->f_graph[it.tp_idx])) {
					value += a->f_Qfunc[it.tp_idx] (thread_id, it.idx, jt.idx, NULL, a->f_Qfunc_arg[it.tp_idx]);
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
					    a->ff_Qfunc[it.tp_idx][jt.tp_idx] (thread_id, it.idx, jt.idx, NULL,
									       a->ff_Qfunc_arg[it.tp_idx][jt.tp_idx]);
				}
			}
			return value;
		}

		case GMRFLib_PREOPT_TP_BETA:
			return value;

		default:
			GMRFLib_ASSERT_RETVAL(0 == 1, GMRFLib_ESNH, 0.0);
		}

		GMRFLib_ASSERT_RETVAL(0 == 1, GMRFLib_ESNH, 0.0);
	}
		break;

	case GMRFLib_PREOPT_TP_BETA:
	{
		switch (jt.tp) {
		case GMRFLib_PREOPT_TP_BETA:
		{
			if (same_tp) {
				value += (a->prior_precision ? a->prior_precision[it.tp_idx] : 0.0);
			}
			return value;
		}

		default:
			GMRFLib_ASSERT_RETVAL(0 == 1, GMRFLib_ESNH, 0.0);
		}

		GMRFLib_ASSERT_RETVAL(0 == 1, GMRFLib_ESNH, 0.0);
	}
		break;

	default:
		GMRFLib_ASSERT_RETVAL(0 == 1, GMRFLib_ESNH, 0.0);
		break;
	}

	return value;
}

double GMRFLib_preopt_like_Qfunc(int thread_id, int node, int nnode, double *UNUSED(values), void *arg)
{
	/*
	 * this is Qfunction for the likelihood part in preopt
	 */
	GMRFLib_preopt_tp *a = (GMRFLib_preopt_tp *) arg;
	GMRFLib_idxval_tp *elm = NULL;
	double *lc = a->like_c[thread_id], value = 0.0;

	if (!lc) {
		return 0.0;
	}
	// imin = node; imax = nnode;
	if (node == nnode) {
		elm = a->AtA_idxval[node][0];
		// value = GMRFLib_dot_product(elm, lc);
		GMRFLib_dot_product_INLINE(value, elm, lc);
	} else {
		int k = 1 + GMRFLib_iwhich_sorted(nnode, a->like_graph->lnbs[node], a->like_graph->lnnbs[node]);
		elm = a->AtA_idxval[node][k];
		// value = GMRFLib_dot_product(elm, lc);
		GMRFLib_dot_product_INLINE(value, elm, lc);
	}

	return value;
}

double GMRFLib_preopt_like_Qfunc_k(int thread_id, int node, int k, double *UNUSED(values), void *arg)
{
	/*
	 * Special version without the need to 'iwhich_sorted' as we know what 'k' is. 
	 */
	GMRFLib_preopt_tp *a = (GMRFLib_preopt_tp *) arg;
	double *lc = a->like_c[thread_id];
	double value = 0.0;

	if (lc) {
		GMRFLib_idxval_tp *elm = a->AtA_idxval[node][k];
		// value = GMRFLib_dot_product(elm, lc);
		GMRFLib_dot_product_INLINE(value, elm, lc);
	}
	return value;
}

double GMRFLib_preopt_Qfunc_OLD(int thread_id, int node, int nnode, double *values, void *arg)
{
	if (nnode < 0) {
		// can be used as a check
		if (!values) {
			return 0.0;
		}
		GMRFLib_preopt_tp *a = (GMRFLib_preopt_tp *) arg;
		int *jj = a->preopt_graph->lnbs[node];
		int lnnbs = a->preopt_graph->lnnbs[node];
		char *like_is_nb = a->preopt_graph_like_is_nb[node];
		char *latent_is_nb = a->preopt_graph_latent_is_nb[node];

		values[0] = GMRFLib_preopt_like_Qfunc(thread_id, node, node, NULL, a->like_Qfunc_arg)
		    + GMRFLib_preopt_latent_Qfunc(thread_id, node, node, NULL, a->latent_Qfunc_arg);
		for (int k = 0; k < lnnbs; k++) {
			values[1 + k] = (like_is_nb[k] ? GMRFLib_preopt_like_Qfunc(thread_id, node, jj[k], NULL, a->like_Qfunc_arg) : 0.0)
			    + (latent_is_nb[k] ? GMRFLib_preopt_latent_Qfunc(thread_id, node, jj[k], NULL, a->latent_Qfunc_arg) : 0.0);
		}
		return 0.0;
	} else {
		GMRFLib_preopt_tp *a = (GMRFLib_preopt_tp *) arg;
		double value = 0.0;

		if (node == nnode) {
			value = GMRFLib_preopt_like_Qfunc(thread_id, node, node, NULL, a->like_Qfunc_arg)
			    + GMRFLib_preopt_latent_Qfunc(thread_id, node, node, NULL, a->latent_Qfunc_arg);
		} else {
			if (GMRFLib_graph_is_nb(node, nnode, a->like_graph)) {
				value = GMRFLib_preopt_like_Qfunc(thread_id, node, nnode, NULL, a->like_Qfunc_arg);
			}
			if (GMRFLib_graph_is_nb(node, nnode, a->latent_graph)) {
				value += GMRFLib_preopt_latent_Qfunc(thread_id, node, nnode, NULL, a->latent_Qfunc_arg);
			}
		}
		return value;
	}
}

double GMRFLib_preopt_Qfunc(int thread_id, int node, int nnode, double *values, void *arg)
{
	if (nnode < 0) {
		// can be used as a check
		if (values) {
			// this used like_Qfunc_k for which we do not need to find the index as we know it already (=kk). we can do this as we
			// do all in a increasing sequence
			GMRFLib_preopt_tp *a = (GMRFLib_preopt_tp *) arg;
			int *jj = a->preopt_graph->lnbs[node];
			int lnnbs = a->preopt_graph->lnnbs[node];
			char *like_is_nb = a->preopt_graph_like_is_nb[node];
			char *latent_is_nb = a->preopt_graph_latent_is_nb[node];

			int kk = 0;
			values[0] = GMRFLib_preopt_like_Qfunc_k(thread_id, node, kk, NULL, a->like_Qfunc_arg)
			    + GMRFLib_preopt_latent_Qfunc(thread_id, node, node, NULL, a->latent_Qfunc_arg);
			for (int k = 0; k < lnnbs; k++) {
				values[1 + k] = (like_is_nb[k] ? GMRFLib_preopt_like_Qfunc_k(thread_id, node, ++kk, NULL, a->like_Qfunc_arg) : 0.0)
				    + (latent_is_nb[k] ? GMRFLib_preopt_latent_Qfunc(thread_id, node, jj[k], NULL, a->latent_Qfunc_arg) : 0.0);
			}
		}
		return 0.0;
	} else {
		GMRFLib_preopt_tp *a = (GMRFLib_preopt_tp *) arg;
		double value = 0.0;

		if (node == nnode) {
			value = GMRFLib_preopt_like_Qfunc(thread_id, node, node, NULL, a->like_Qfunc_arg)
			    + GMRFLib_preopt_latent_Qfunc(thread_id, node, node, NULL, a->latent_Qfunc_arg);
		} else {
			if (GMRFLib_graph_is_nb(node, nnode, a->like_graph)) {
				value = GMRFLib_preopt_like_Qfunc(thread_id, node, nnode, NULL, a->like_Qfunc_arg);
			}
			if (GMRFLib_graph_is_nb(node, nnode, a->latent_graph)) {
				value += GMRFLib_preopt_latent_Qfunc(thread_id, node, nnode, NULL, a->latent_Qfunc_arg);
			}
		}
		return value;
	}
}

double GMRFLib_preopt_gcpo_Qfunc(int thread_id, int node, int nnode, double *UNUSED(values), void *arg)
{
	// this function is special. the graph is preopt, but only the prior is returned.

	if (nnode < 0) {
		return NAN;
	}

	GMRFLib_preopt_tp *a = (GMRFLib_preopt_tp *) arg;
	double value = 0.0;

	if (node == nnode || GMRFLib_graph_is_nb(node, nnode, a->latent_graph)) {
		value = GMRFLib_preopt_latent_Qfunc(thread_id, node, nnode, NULL, a->latent_Qfunc_arg);
		if (a->gcpo_mask) {
			value *= a->gcpo_mask[node] * a->gcpo_mask[nnode];
			if (node == nnode) {
				value += a->gcpo_diag[node];
			}
		}
	}

	return value;
}

double GMRFLib_preopt_Qfunc_like(int thread_id, int node, int nnode, double *UNUSED(values), void *arg)
{
	// standalone function to return the likelihood part only
	if (nnode < 0) {
		return NAN;
	}

	GMRFLib_preopt_tp *a = (GMRFLib_preopt_tp *) arg;
	double value = 0.0;
	if (node == nnode || GMRFLib_graph_is_nb(node, nnode, a->like_graph)) {
		value = GMRFLib_preopt_like_Qfunc(thread_id, node, nnode, NULL, a->like_Qfunc_arg);
	}

	return value;
}

double GMRFLib_preopt_Qfunc_prior(int thread_id, int node, int nnode, double *UNUSED(values), void *arg)
{
	// standalone function to return the prior part
	if (nnode < 0) {
		return NAN;
	}

	GMRFLib_preopt_tp *a = (GMRFLib_preopt_tp *) arg;
	double value = 0.0;
	if (node == nnode || GMRFLib_graph_is_nb(node, nnode, a->latent_graph)) {
		value = GMRFLib_preopt_latent_Qfunc(thread_id, node, nnode, NULL, a->latent_Qfunc_arg);
	}

	return value;
}


int GMRFLib_preopt_bnew(int thread_id, double *b, GMRFLib_preopt_tp *preopt)
{
	GMRFLib_ENTER_ROUTINE;
	GMRFLib_preopt_bnew_like(b, preopt->like_b[thread_id], preopt);
	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}

int GMRFLib_preopt_bnew_like(double *bnew, double *blike, GMRFLib_preopt_tp *preopt)
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

	// bnew[i] = GMRFLib_dot_product(elm, blike);

#define CODE_BLOCK							\
	for (int i = 0; i < preopt->n; i++) {				\
		if (A[i]) {						\
			GMRFLib_idxval_tp *elm = A[i];			\
			GMRFLib_dot_product_INLINE(bnew[i], elm, blike); \
		}							\
	}

	RUN_CODE_BLOCK(GMRFLib_MAX_THREADS(), 0, 0);
#undef CODE_BLOCK

	return GMRFLib_SUCCESS;
}

int GMRFLib_preopt_predictor(double *predictor, double *latent, GMRFLib_preopt_tp *preopt)
{
	GMRFLib_ENTER_ROUTINE;
	int val = GMRFLib_preopt_predictor_core(predictor, latent, preopt, 1);
	GMRFLib_LEAVE_ROUTINE;
	return val;
}

int GMRFLib_preopt_full_predictor(double *predictor, double *latent, GMRFLib_preopt_tp *preopt)
{
	GMRFLib_ENTER_ROUTINE;
	int val = GMRFLib_preopt_predictor_core(predictor, latent, preopt, 0);
	GMRFLib_LEAVE_ROUTINE;
	return val;
}

int GMRFLib_preopt_predictor_core(double *predictor, double *latent, GMRFLib_preopt_tp *preopt, int likelihood_only)
{
	// if likelihood_only, only compute the part that is needed for the likelihood.

	// if !likelihood_only, compute the whole predictor

	double *pred = Calloc(preopt->mnpred, double);
	int data_rich_case = GMRFLib_preopt_predictor_strategy;
	// int data_rich_case = (IMAX(preopt->mpred, preopt->npred) > preopt->n);
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
			double *pred_offset = pred + offset;
#define CODE_BLOCK							\
			for (int j = 0; j < 2; j++) {			\
				if (j == 0) {				\
					for (int i = 0; i < preopt->n; i++) { \
						GMRFLib_idxval_tp *At = preopt->At_idxval[i]; \
						double lat = latent[i];	\
						for (int k = 0; k < At->n; k++) { \
							pred_offset[At->idx[k]] += lat * At->val[k]; \
						}			\
					}				\
				} else {				\
					for (int i = 0; i < preopt->n; i++) { \
						GMRFLib_idxval_tp *pAAt = preopt->pAAt_idxval[i]; \
						double lat = latent[i];	\
						for (int k = 0; k < pAAt->n; k++) { \
							pred[pAAt->idx[k]] += lat * pAAt->val[k]; \
						}			\
					}				\
				}					\
			}

			RUN_CODE_BLOCK(2, 0, 0);
#undef CODE_BLOCK
		} else {
			// one loop
			double *pred_offset = pred + offset;
			for (int i = 0; i < preopt->n; i++) {
				GMRFLib_idxval_tp *At = preopt->At_idxval[i];
				double lat = latent[i];
				for (int k = 0; k < At->n; k++) {
					pred_offset[At->idx[k]] += lat * At->val[k];
				}
			}
		}

	} else {

		// not data-rich case

		if (preopt->pA_idxval) {
			// both loops
			double *pred_offset = pred + offset;

			// pred_offset[i] = GMRFLib_dot_product(elm, latent); 
			// pred[i] = GMRFLib_dot_product(elm, latent); 

#define CODE_BLOCK							\
			for (int j = 0; j < 2; j++) {			\
				if (j == 0) {				\
					for (int i = 0; i < preopt->npred; i++) { \
						GMRFLib_idxval_tp *elm = preopt->A_idxval[i]; \
						GMRFLib_dot_product_INLINE(pred_offset[i], elm, latent); \
					}				\
				} else {				\
					for (int i = 0; i < preopt->mpred; i++) { \
						GMRFLib_idxval_tp *elm = preopt->pAA_idxval[i]; \
						GMRFLib_dot_product_INLINE(pred[i], elm, latent); \
					}				\
				}					\
			}

			RUN_CODE_BLOCK(2, 0, 0);
#undef CODE_BLOCK

		} else {
			// one loop
			double *pred_offset = pred + offset;

			// pred_offset[i] = GMRFLib_dot_product(elm, latent); 
#define CODE_BLOCK							\
			for (int i = 0; i < preopt->npred; i++) {	\
				if (preopt->A_idxval[i]) {		\
					GMRFLib_idxval_tp *elm = preopt->A_idxval[i]; \
					GMRFLib_dot_product_INLINE(pred_offset[i], elm, latent); \
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

	return GMRFLib_SUCCESS;
}

int GMRFLib_preopt_predictor_moments(double *mean, double *variance, GMRFLib_preopt_tp *preopt, GMRFLib_problem_tp *problem, double *optional_mean)
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
				// mean[i] = GMRFLib_dot_product(elm, mm); 
#define CODE_BLOCK							\
				for (int i = 0; i < mpred; i++) {	\
					GMRFLib_idxval_tp *elm = preopt->pAA_idxval[i]; \
					GMRFLib_dot_product_INLINE(mean[i], elm, mm); \
				}

				RUN_CODE_BLOCK(GMRFLib_MAX_THREADS(), 0, 0);
#undef CODE_BLOCK
			}
		} else {
			// both mean and variance

#define CODE_BLOCK							\
			for (int i = 0; i < mpred; i++) {		\
				double m = 0.0, var = 0.0, *cov = NULL;	\
				int k, j, kk, jj;			\
				GMRFLib_idxval_tp *elm = preopt->pAA_idxval[i]; \
				for (k = 0; k < preopt->pAA_idxval[i]->n; k++) {	\
					j = elm->idx[k];		\
					if (compute_mean) {		\
						m += elm->val[k] * mm[j]; \
					}				\
					cov = GMRFLib_Qinv_get(problem, j, j); \
					var += SQR(elm->val[k]) * *cov;	\
					double tvar = 0.0;		\
					for (kk = k+1; kk < preopt->pAA_idxval[i]->n; kk++){ \
						jj = elm->idx[kk];	\
						cov = GMRFLib_Qinv_get(problem, j, jj);	\
						tvar += elm->val[kk] * *cov; \
					}				\
					var += 2.0 * elm->val[k] * tvar; \
				}					\
				if (compute_mean) {			\
					mean[i] = m;			\
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
		double *mean_offset = mean + offset;
		if (data_rich_case) {
			for (int i = 0; i < preopt->n; i++) {
				GMRFLib_idxval_tp *At = preopt->At_idxval[i];
				double lat = mm[i];
				for (int k = 0; k < At->n; k++) {
					mean_offset[At->idx[k]] += lat * At->val[k];
				}
			}
		} else {
			// mean_offset[i] += GMRFLib_dot_product(elm, mm); 
#define CODE_BLOCK							\
			for (int i = 0; i < npred; i++) {		\
				GMRFLib_idxval_tp *elm = preopt->A_idxval[i]; \
				GMRFLib_dot_product_INLINE_ADDTO(mean_offset[i], elm, mm); \
			}

			RUN_CODE_BLOCK(GMRFLib_MAX_THREADS(), 0, 0);
#undef CODE_BLOCK
		}
	} else {

		// both mean and variance
		double *mean_offset = (compute_mean ? mean + offset : NULL);
		double *variance_offset = variance + offset;

#define CODE_BLOCK							\
		for (int i = 0; i < npred; i++) {			\
			double m = 0.0, var = 0.0, zero = 0.0, *cov = NULL; \
			int k, j, kk, jj;				\
			GMRFLib_idxval_tp *elm = preopt->A_idxval[i];	\
			for (k = 0; k < preopt->A_idxval[i]->n; k++){	\
				j = elm->idx[k];			\
				if (compute_mean) {			\
					m += elm->val[k] * mm[j];	\
				}					\
				cov = GMRFLib_Qinv_get(problem, j, j);	\
				var += SQR(elm->val[k]) * *cov;		\
				double tvar = 0.0;			\
				for (kk = k+1; kk < preopt->A_idxval[i]->n; kk++){ \
					jj = elm->idx[kk];		\
					cov = GMRFLib_Qinv_get(problem, j, jj);	\
					if (!cov) {			\
						err_count++; /* its ok not todo critical */ \
						cov = &zero;		\
					}				\
					tvar += elm->val[kk] * *cov;	\
				}					\
				var += 2.0 * elm->val[k] * tvar;	\
			}						\
			if (compute_mean) {				\
				mean_offset[i] = m;			\
			}						\
			variance_offset[i] = var;			\
		}

		RUN_CODE_BLOCK(GMRFLib_MAX_THREADS(), 0, 0);
#undef CODE_BLOCK
	}

	if (err_count) {
		static int shown = 0;
		if (!shown) {
			shown = 1;
			fprintf(stdout, "\n\n%s:%d:(%s)\n\tMissing (about) %1d covariances.\n\t%s\n\t%s\n\t%s\n\n\n",
				__FILE__, __LINE__, __GMRFLib_FuncName, err_count,
				"Either the A-matrix has not the required rank",
				"or correlations are numerically zero.", "Further warnings are disabled.");
		}
	}

	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}

int GMRFLib_preopt_update(int thread_id, GMRFLib_preopt_tp *preopt, double *like_b, double *like_c)
{
	GMRFLib_ENTER_ROUTINE;
	int np = preopt->Npred;

	if (!(preopt->like_b[thread_id])) {
		preopt->like_b[thread_id] = Calloc(np, double);
		preopt->like_c[thread_id] = Calloc(np, double);
	}
	Memcpy(preopt->like_b[thread_id], like_b, np * sizeof(double));
	Memcpy(preopt->like_c[thread_id], like_c, np * sizeof(double));

	if (!(preopt->total_b[thread_id])) {
		preopt->total_b[thread_id] = Calloc(preopt->n, double);
	}
	Memset(preopt->total_b[thread_id], 0, preopt->n * sizeof(double));
	GMRFLib_preopt_bnew(thread_id, preopt->total_b[thread_id], preopt);

	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}

int GMRFLib_preopt_free(GMRFLib_preopt_tp *preopt)
{
	if (!preopt) {
		return GMRFLib_SUCCESS;
	}
#pragma omp parallel sections
	{
#pragma omp section
		{
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

			for (int i = 0; i < preopt->preopt_graph->n; i++) {
				Free(preopt->preopt_graph_latent_is_nb[i]);
				Free(preopt->preopt_graph_like_is_nb[i]);
			}
			Free(preopt->preopt_graph_latent_is_nb);
			Free(preopt->preopt_graph_like_is_nb);

			GMRFLib_graph_free(preopt->preopt_graph);
			GMRFLib_graph_free(preopt->like_graph);
			GMRFLib_graph_free(preopt->latent_graph);
			GMRFLib_free_constr(preopt->latent_constr);
		}
	}
	Free(preopt);

	return GMRFLib_SUCCESS;
}

double *GMRFLib_preopt_measure_time(int thread_id, GMRFLib_preopt_tp *preopt, double *res, double *test_vector)
{
	// return alloc'ed double *cpu measurements.
	// cpu[0] and cpu[1] is the time for doing Q %*% x.
	// cpu[0] does elementwise computations, cpu[1] does blockwise using tabulated Q.

	double value = 0.0;
	double *cpu = Calloc(3, double);

	GMRFLib_graph_tp *like_graph = preopt->like_graph;
	void *like_Qfunc_arg = preopt->like_Qfunc_arg;

	// this will be measure with serial or with group
	cpu[0] = -GMRFLib_timer();
	for (int i = 0; i < like_graph->n; i++) {
		value += GMRFLib_preopt_like_Qfunc(thread_id, i, i, NULL, like_Qfunc_arg);
		for (int jj = 0, j; jj < like_graph->lnnbs[i]; jj++) {
			j = like_graph->lnbs[i][jj];
			value += GMRFLib_preopt_like_Qfunc(thread_id, i, j, NULL, like_Qfunc_arg);
		}
	}
	cpu[0] += GMRFLib_timer();
	assert(!ISNAN(value));

	GMRFLib_Qfunc_tp *Qfunc = preopt->preopt_Qfunc;
	GMRFLib_graph_tp *graph = preopt->preopt_graph;
	void *Qfunc_arg = preopt->preopt_Qfunc_arg;

	Calloc_init(2 * graph->n, 2);
	double *x = NULL;
	double *xx = Calloc_get(graph->n);
	if (!test_vector) {
		x = Calloc_get(graph->n);
		for (int i = 0; i < graph->n; i++) {
			x[i] = GMRFLib_uniform();
		}
	} else {
		x = test_vector;
	}

	GMRFLib_tabulate_Qfunc_tp *tab = NULL;
	GMRFLib_tabulate_Qfunc_core(thread_id, &tab, graph, Qfunc, Qfunc_arg, NULL, 1);

	// this will be measured with serial or parallel
	cpu[1] = -GMRFLib_timer();
	GMRFLib_Qx(thread_id, xx, x, graph, tab->Qfunc, tab->Qfunc_arg);
	cpu[1] += GMRFLib_timer();

	if (0) {
		double check = 0.0;
		for (int i = 0; i < graph->n; i++) {
			check += ABS(xx[i]);
		}
	}

	if (res) {
		int inc = 1;
		res[0] = value;
		res[1] = dasum_(&(graph->n), xx, &inc);
	}

	Calloc_free();
	GMRFLib_free_tabulate_Qfunc(tab);

	return cpu;
}

double *GMRFLib_preopt_measure_time2(GMRFLib_preopt_tp *preopt)
{
	// cpu[0] measure predictor calculations for !data_rich and data_rich case
	double *cpu = Calloc(1, double);

	Calloc_init(preopt->mnpred + preopt->n, 2);
	double *pred = Calloc_get(preopt->mnpred);
	double *lat = Calloc_get(preopt->n);

	for (int i = 0; i < preopt->n; i++) {
		lat[i] = GMRFLib_uniform();
	}
	cpu[0] = -GMRFLib_timer();
	GMRFLib_preopt_full_predictor(pred, lat, preopt);
	cpu[0] += GMRFLib_timer();
	Calloc_free();

	return cpu;
}


/* hgmrfm.c
 * 
 * Copyright (C) 2007-2024 Havard Rue
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

#include <time.h>
#include <strings.h>
#include <stdlib.h>
#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

#define HYPER_NEW(name_, initial_)					\
	if (1) {							\
		name_ = Calloc(GMRFLib_MAX_THREADS(), double *);	\
		for(int i_ = 0; i_ < GMRFLib_MAX_THREADS(); i_++) {	\
			name_[i_] = Calloc(1, double);			\
			name_[i_][0] = initial_;			\
		}							\
	}

int GMRFLib_init_hgmrfm(GMRFLib_hgmrfm_tp **hgmrfm, int n, int n_ext,
			int *eta_sumzero, double **logprec_unstruct_omp,
			const char *Aext_fnm, double Aext_precision,
			int nf, int **c, double **w,
			GMRFLib_graph_tp **f_graph, GMRFLib_Qfunc_tp **f_Qfunc,
			void **f_Qfunc_arg, char *f_sumzero, GMRFLib_constr_tp **f_constr,
			GMRFLib_Qfunc_tp ***ff_Qfunc, void ***ff_Qfunc_arg,
			int nbeta, double **covariate, double *prior_precision, GMRFLib_ai_param_tp *UNUSED(ai_par))
{
	/*
	 * define a HGMRF-model, of the form
	 * 
	 * \y_i ~ f(eta_i)
	 * 
	 * and
	 * 
	 * \eta_i | ... ~ N( \sum_k f_kc_k(i) + \sum_m z_mi\beta_m, 1/\lambda_unstruct)
	 *
	 *
	 * About the interaction terms: ThHe ff_Qfunc is a 2d array [0...nf-1] x [0...nf-1] of interaction functions, which adds interaction terms between f[i][k]
	 * and f[j][k], for i!=j and k = 0...MIN(f_graph[i]->n, f_graph[j]->n). Yes, it's a bit weird, but with this feature we can circumwent the requirement that
	 * one field is only allowed to be present once in the predictor. This feature is for internal use only, really...
	 */

#define SET_ELEMENT(i_, j_, Qij_, id_) SET_ELEMENT_ADV(i_, j_, Qij_, 0, id_)
#define SET_ELEMENT_FORCE(i_, j_, Qij_, id_) SET_ELEMENT_ADV(i_, j_, Qij_, 1, id_)
#define SET_ELEMENT_ADV(i_, j_, Qij_, test, id_) \
	if (Qij_ || (test)) {						\
		if (ntriples[id_] >= ntriples_max[id_]) {		\
			ntriples_max[id_] += n;				\
			ilist[id_] = Realloc(ilist[id_], ntriples_max[id_], int); \
			jlist[id_] = Realloc(jlist[id_], ntriples_max[id_], int); \
			Qijlist[id_] = Realloc(Qijlist[id_], ntriples_max[id_], double); \
		}							\
		ilist[id_][ntriples[id_]] = i_;				\
		jlist[id_][ntriples[id_]] = j_;				\
		Qijlist[id_][ntriples[id_]] = Qij_;			\
		ntriples[id_]++;					\
	}

#define SET_ELEMENT_LC(i_, j_, Qij_) SET_ELEMENT_ADV_LC(i_, j_, Qij_, 0)
#define SET_ELEMENT_FORCE_LC(i_, j_, Qij_) SET_ELEMENT_ADV_LC(i_, j_, Qij_, 1)
#define SET_ELEMENT_ADV_LC(i_, j_, Qij_, test) if (Qij_ || (test)) {	\
		if (ntriples_lc >= ntriples_max_lc) {			\
			ntriples_max_lc += n;				\
			ilist_lc = Realloc(ilist_lc, ntriples_max_lc, int); \
			jlist_lc = Realloc(jlist_lc, ntriples_max_lc, int); \
			Qijlist_lc = Realloc(Qijlist_lc, ntriples_max_lc, double); \
		}							\
		if(0)printf("set_element_lc %d %d %g\n", i_, j_, Qij_);	\
		ilist_lc[ntriples_lc] = i_;				\
		jlist_lc[ntriples_lc] = j_;				\
		Qijlist_lc[ntriples_lc] = Qij_;				\
		ntriples_lc++;						\
	}

	int i, ii, j, jj, k, l, m, N, **ilist = NULL, **jlist = NULL, *ntriples = NULL, *ntriples_max = NULL,
	    idx_map_eta = 0, *idx_map_f = NULL, *idx_map_beta = NULL, *idx_map_lc = NULL, offset, ***fidx = NULL, **nfidx = NULL, **lfidx =
	    NULL, fidx_add = 5;
	int nu = 0, *uniq = NULL;
	double **Qijlist = NULL, value, **ww = NULL;
	GMRFLib_hgmrfm_arg_tp *arg = NULL;
	GMRFLib_constr_tp *fc = NULL;

	if (!hgmrfm) {
		return GMRFLib_SUCCESS;
	}

	*hgmrfm = Calloc(1, GMRFLib_hgmrfm_tp);
	arg = Calloc(1, GMRFLib_hgmrfm_arg_tp);
	n = IMAX(0, n);
	nbeta = IMAX(0, nbeta);
	nf = IMAX(0, nf);
	arg->n = n;
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
						GMRFLib_ASSERT(ff_Qfunc_arg[i][j] == ff_Qfunc_arg[j][i], GMRFLib_EPARAMETER);
					}
				}
			}
		}
	}

	if (nf) {
		/*
		 * make the weights for the computation, taken into account that the default weight is 1.0 
		 */
		ww = Calloc(nf, double *);
		for (k = 0; k < nf; k++) {
			ww[k] = Calloc(n, double);
			if (w && w[k]) {
				Memcpy(ww[k], w[k], n * sizeof(double));
			} else {
				for (j = 0; j < n; j++) {
					ww[k][j] = 1.0;
				}
			}
		}
	}

	/*
	 * need to read the Aext matrix, which is required. Allocate a new copy of precision, as we need to make it persistent. 
	 */
	if (Aext_fnm) {
		double **lprec_omp;
		HYPER_NEW(lprec_omp, log(Aext_precision));
		GMRFLib_tabulate_Qfunc_from_file(&(arg->eta_ext_Q), &(arg->eta_ext_graph), Aext_fnm, -1, lprec_omp);
		GMRFLib_ASSERT(arg->eta_ext_graph->n == n + n_ext, GMRFLib_EPARAMETER);	/* this is required!!!!! */
		arg->n_ext = n_ext;
		// GMRFLib_printf_graph(stdout, arg->eta_ext_graph);
	} else {
		arg->eta_ext_Q = NULL;
		arg->eta_ext_graph = NULL;
		arg->n_ext = 0;
	}

	/*
	 * Our first job, is to go through the model and compute all interactions etc that are defined through the \eta-model. 
	 * define the index-mapping. The outline of x is (eta, f[0], ..., beta[0], ...)
	 */
	idx_map_eta = offset = arg->n_ext;		       /* starting for 'eta' (not eta_ext) */
	offset += n;					       /* including eta */

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

	/*
	 * we need to make sure that all nodes are present in the eta-graph. We do this by adding just zero's. This is not checked for when
	 * building the graph, but we do it here, so we use 1 as last argument to SET_ELEMENT(,,1). 
	 */

	int tmax = GMRFLib_MAX_THREADS();

	ilist = Calloc(tmax, int *);
	jlist = Calloc(tmax, int *);
	Qijlist = Calloc(tmax, double *);
	ntriples = Calloc(tmax, int);
	ntriples_max = Calloc(tmax, int);

	/*
	 * set for all threads 
	 */
	N = offset;					       /* N is the grand-total. */
	SET_ELEMENT_FORCE(N - 1, N - 1, 0.0, 0);

	/*
	 * Ensure also that eta_ext elements are set to zero, so we produce the joint graph. yes, we need the hole graph!!!! 
	 */
	if (idx_map_eta > 0) {
		assert(arg->eta_ext_graph);

		GMRFLib_graph_tp *g = arg->eta_ext_graph;
		for (i = 0; i < g->n; i++) {
			SET_ELEMENT_FORCE(i, i, 0.0, 0);
			for (j = 0; j < g->nnbs[i]; j++) {
				jj = g->nbs[i][j];
				SET_ELEMENT_FORCE(i, jj, 0.0, 0);
			}
		}
	}

	/*
	 * If we have cross-terms, make sure to mark these cross-terms as neigbours; just fill them with zero's. If they are connected in the data, then this value
	 * will be overrided and this is how it should be.
	 */
	if (ff_Qfunc) {
		for (i = 0; i < nf; i++) {
			for (j = i + 1; j < nf; j++) {
				if (ff_Qfunc[i][j]) {
					for (k = 0; k < IMIN(f_graph[i]->n, f_graph[j]->n); k++) {
						SET_ELEMENT_FORCE(idx_map_f[i] + k, idx_map_f[j] + k, 0.0, 0);
					}
				}
			}
		}
	}

	if (nf) {
		/*
		 * this is required for making the ffield computations fast. we need to build an index array for the ffields:
		 * fidx[k][m] is a list of length nidx[k][m] containing all those i's where c[k][i]=m. lfidx[k][m] is just the
		 * alloced length that is grown dynamically. 
		 */
		fidx = Calloc(nf, int **);
		nfidx = Calloc(nf, int *);
		lfidx = Calloc(nf, int *);
		for (k = 0; k < nf; k++) {
			fidx[k] = Calloc(f_graph[k]->n, int *);
			nfidx[k] = Calloc(f_graph[k]->n, int);
			lfidx[k] = Calloc(f_graph[k]->n, int);
			for (m = 0; m < f_graph[k]->n; m++) {
				nfidx[k][m] = 0;
				lfidx[k][m] = 0;	       /* initalise to zero length to minimise storage */
				fidx[k][m] = NULL;
			}
		}
		for (k = 0; k < nf; k++) {
			int lenf = f_graph[k]->n;
			for (i = 0; i < n; i++) {
				m = c[k][i];
				if (LEGAL(m, lenf)) {
					if (nfidx[k][m] >= lfidx[k][m]) {
						lfidx[k][m] += fidx_add;
						fidx[k][m] = Realloc(fidx[k][m], lfidx[k][m], int);
					}
					fidx[k][m][nfidx[k][m]] = i;
					nfidx[k][m]++;
				}
			}
		}
	}

	int max_threads = IMIN(5, GMRFLib_openmp->max_threads_outer);
#pragma omp parallel sections num_threads(max_threads)
	{
#pragma omp section
		{
			/*
			 * \eta_i^2 = 1 
			 */
			int thread = omp_get_thread_num();
			int it;

			for (it = 0; it < n; it++) {
				SET_ELEMENT(idx_map_eta + it, idx_map_eta + it, 1.0, thread);
			}
		}
#pragma omp section
		{
			/*
			 * \eta_i f_jk = - 1_{c_j(i) = k} 
			 */
			if (nf) {
				int thread = omp_get_thread_num();
				int jt, kt, iit, it;

				for (jt = 0; jt < nf; jt++) {
					for (kt = 0; kt < f_graph[jt]->n; kt++) {
						for (iit = 0; iit < nfidx[jt][kt]; iit++) {
							it = fidx[jt][kt][iit];
							SET_ELEMENT(idx_map_eta + it, idx_map_f[jt] + kt, -ww[jt][it], thread);
						}
					}
				}
			}
		}

#pragma omp section
		{
			/*
			 * \eta_i \beta_j 
			 */
			if (nbeta) {
				int thread = omp_get_thread_num();
				int jt, it;

				for (jt = 0; jt < nbeta; jt++) {
					for (it = 0; it < n; it++) {
						SET_ELEMENT(idx_map_eta + it, idx_map_beta[jt], -covariate[jt][it], thread);
					}
				}
			}
		}
#pragma omp section
		{
			/*
			 * f_jk beta_m = \sum z_ki, for all i: c_j[i] = k 
			 */
			if (nf && nbeta) {
				int thread = omp_get_thread_num();
				int jt, kt, mt, it, iit;

				for (jt = 0; jt < nf; jt++) {
					for (kt = 0; kt < f_graph[jt]->n; kt++) {
						for (mt = 0; mt < nbeta; mt++) {
							double valuet = 0.0;

							for (iit = 0; iit < nfidx[jt][kt]; iit++) {
								it = fidx[jt][kt][iit];
								valuet += covariate[mt][it] * ww[jt][it];
							}
							SET_ELEMENT(idx_map_f[jt] + kt, idx_map_beta[mt], valuet, thread);
						}
					}
				}
			}
		}

#pragma omp section
		{
			/*
			 * beta_k beta_m = sum_i z_ki z_mi 
			 */
			if (nbeta) {
				int thread = omp_get_thread_num();
				int kt, mt, it;

				for (kt = 0; kt < nbeta; kt++) {
					for (mt = kt; mt < nbeta; mt++) {
						double valuet = 0.0;

						for (it = 0; it < n; it++) {
							valuet += covariate[kt][it] * covariate[mt][it];
						}
						SET_ELEMENT(idx_map_beta[kt], idx_map_beta[mt], valuet, thread);
					}
				}
			}
		}
	}						       /* END of parallel sections */

	/*
	 * this one could be heavy... 
	 */

	/*
	 * f_jk f_ml = \sum 1_{i : c_j(i) == k && c_m(i) == l } 
	 */
	if (nf) {
		/*
		 * first we build a 1d index to emulate a 2d index. 
		 */
		int jm, jm_idx;
		int *j_idx = Calloc(ISQR(nf), int);
		int *m_idx = Calloc(ISQR(nf), int);

		for (j = 0, jm = 0; j < nf; j++) {
			for (m = j; m < nf; m++) {
				j_idx[jm] = j;
				m_idx[jm] = m;
				jm++;
			}
		}

		/*
		 * try go get better load-balancing moving the parallell into an inner loop (NEW) instead of making the first loop parallel (OLD) 
		 */
		if (1) {
			/*
			 * NEW VERSION, parallel inner loop 
			 */

			// FIXME("************************NEW");

			for (jm_idx = 0; jm_idx < jm; jm_idx++) {
				j = j_idx[jm_idx];
				m = m_idx[jm_idx];

#pragma omp parallel for private(k, l, value, ii, i) num_threads(GMRFLib_openmp->max_threads_outer)
				for (k = 0; k < f_graph[j]->n; k++) {
					int thread = omp_get_thread_num();
					for (l = 0; l < f_graph[m]->n; l++) {
						value = 0.0;
						if (nfidx[j][k] < nfidx[m][l]) {
							for (ii = 0; ii < nfidx[j][k]; ii++) {
								i = fidx[j][k][ii];
								if (c[m][i] == l) {
									value += ww[j][i] * ww[m][i];
								}
							}
						} else {
							for (ii = 0; ii < nfidx[m][l]; ii++) {
								i = fidx[m][l][ii];
								if (c[j][i] == k) {
									value += ww[j][i] * ww[m][i];
								}
							}
						}
						SET_ELEMENT(idx_map_f[j] + k, idx_map_f[m] + l, value, thread);
					}
				}
			}
		} else {
			/*
			 * OLD VERSION, parallel outher loop 
			 */

			// FIXME("*****************OLD");

#pragma omp parallel for private(jm_idx, j, k, m, l, value, ii, i) num_threads(GMRFLib_openmp->max_threads_outer)
			for (jm_idx = 0; jm_idx < jm; jm_idx++) {
				int thread = omp_get_thread_num();

				j = j_idx[jm_idx];
				m = m_idx[jm_idx];

				for (k = 0; k < f_graph[j]->n; k++) {
					for (l = 0; l < f_graph[m]->n; l++) {
						value = 0.0;
						if (nfidx[j][k] < nfidx[m][l]) {
							for (ii = 0; ii < nfidx[j][k]; ii++) {
								i = fidx[j][k][ii];
								if (c[m][i] == l) {
									value += ww[j][i] * ww[m][i];
								}
							}
						} else {
							for (ii = 0; ii < nfidx[m][l]; ii++) {
								i = fidx[m][l][ii];
								if (c[j][i] == k) {
									value += ww[j][i] * ww[m][i];
								}
							}
						}
						SET_ELEMENT(idx_map_f[j] + k, idx_map_f[m] + l, value, thread);
					}
				}
			}
		}
	}

	if (0) {
		for (k = 0; k < tmax; k++)
			for (i = 0; i < ntriples[k]; i++) {
				printf("QQ %d %d %d %g\n", k, ilist[k][i], jlist[k][i], Qijlist[k][i]);
			}
	}

	/*
	 * now we need to collect all results... 
	 */
	int *iilist = NULL;
	int *jjlist = NULL;
	double *QQijlist = NULL;
	int nntriples;

	nntriples = 0;
	for (i = 0; i < tmax; i++) {
		nntriples += ntriples[i];
	}

	if (0) {
		for (i = 0; i < tmax; i++) {
			printf("i %d ntriples %d\n", i, ntriples[i]);
		}
	}

	iilist = Calloc(nntriples, int);
	jjlist = Calloc(nntriples, int);
	QQijlist = Calloc(nntriples, double);

	for (i = 0, k = 0; i < tmax; i++) {
		int len = ntriples[i];

		if (len) {
			Memcpy(&iilist[k], ilist[i], len * sizeof(int));
			Memcpy(&jjlist[k], jlist[i], len * sizeof(int));
			Memcpy(&QQijlist[k], Qijlist[i], len * sizeof(double));
			k += len;
		}
	}

	if (0) {
		for (i = 0; i < nntriples; i++) {
			printf("ALL %d %d %g\n", iilist[i], jjlist[i], QQijlist[i]);
		}
	}
	GMRFLib_tabulate_Qfunc_from_list(&(arg->eta_Q), &(arg->eta_graph), nntriples, iilist, jjlist, QQijlist, -1, logprec_unstruct_omp);
	/*
	 * cleanup already here, as we do not need them anymore 
	 */
	Free(iilist);
	Free(jjlist);
	Free(QQijlist);
	for (i = 0; i < tmax; i++) {
		Free(ilist[i]);
		Free(jlist[i]);
		Free(Qijlist[i]);
	}
	Free(ilist);
	Free(jlist);
	Free(Qijlist);
	Free(ntriples);
	Free(ntriples_max);

	/*
	 * now it is time to create the full graph by inserting the ones from the ffields.
	 */
	GMRFLib_ged_tp *ged = NULL;
	GMRFLib_ged_init(&ged, arg->eta_graph);
	if (nf) {
		for (j = 0; j < nf; j++) {
			GMRFLib_ged_insert_graph(ged, f_graph[j], idx_map_f[j]);
		}
	}
	GMRFLib_ged_build(&((*hgmrfm)->graph), ged);
	GMRFLib_ged_free(ged);
	(*hgmrfm)->Qfunc = GMRFLib_hgmrfm_Qfunc;
	(*hgmrfm)->Qfunc_arg = (void *) arg;

	/*
	 * build the constraint, if any. Only simple sum-to-zero constraints are supported.
	 */
	int nconstr = 0;
	if (nf && f_sumzero) {
		for (k = 0; k < nf; k++) {
			nconstr += (f_sumzero[k] ? 1 : 0);
		}
	}

	GMRFLib_iuniques(&nu, &uniq, eta_sumzero, n + n_ext);
	nconstr += nu;

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
		constr->a_matrix = Calloc((*hgmrfm)->graph->n * nconstr, double);
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

		if (nu) {
			for (k = 0; k < nu; k++) {
				ii = uniq[k];
				for (i = 0; i < n + n_ext; i++) {
					constr->a_matrix[i * constr->nc + constr_no] = (eta_sumzero[i] == ii ? 1.0 : 0.0);
				}
				constr->e_vector[constr_no] = 0.0;
				constr_no++;
			}
		}

		GMRFLib_ASSERT(constr_no == constr->nc, GMRFLib_ESNH);
		GMRFLib_prepare_constr(constr, (*hgmrfm)->graph, 0);
		(*hgmrfm)->constr = constr;
	} else {
		(*hgmrfm)->constr = NULL;
	}
	arg->idx_map_f = idx_map_f;
	arg->idx_map_beta = idx_map_beta;
	arg->idx_map_lc = idx_map_lc;
	arg->N = (*hgmrfm)->graph->n;
	GMRFLib_ASSERT(arg->N == N, GMRFLib_ESNH);
	arg->what_type = Calloc(N, GMRFLib_hgmrfm_type_tp);
	for (i = 0; i < N; i++) {
		arg->what_type[i] = GMRFLib_hgmrfm_what_type(i, arg);
	}

	Free(uniq);
	if (nf) {
		for (k = 0; k < nf; k++) {
			for (m = 0; m < f_graph[k]->n; m++) {
				Free(fidx[k][m]);
			}
			Free(fidx[k]);
			Free(nfidx[k]);
			Free(lfidx[k]);
		}
		Free(fidx);
		Free(nfidx);
		Free(lfidx);
		for (k = 0; k < nf; k++) {
			Free(ww[k]);
		}
		Free(ww);
	}
#undef SET_ELEMENT
#undef SET_ELEMENT_LC
#undef SET_ELEMENT_FORCE
#undef SET_ELEMENT_FORCE_LC
#undef SET_ELEMENT_ADV
#undef SET_ELEMENT_ADV_LC

	if (0) {
		GMRFLib_hgmrfm_tp *h = *hgmrfm;

		printf("view hgmrf\n");
		GMRFLib_printf_graph(stdout, h->graph);

		int nn = h->graph->n;
		if (h->constr && h->constr->nc) {
			for (j = 0; j < h->constr->nc; j++) {
				printf("constr %d\n", j);

				for (i = 0; i < nn; i++) {
					printf("%.1f ", h->constr->a_matrix[i * h->constr->nc + j]);
				}
				printf("| %.f\n", h->constr->e_vector[j]);
			}
		}

		int thread_id = 0;
		assert(omp_get_thread_num() == 0);
		GMRFLib_printf_Qfunc(thread_id, stdout, h->graph, h->Qfunc, h->Qfunc_arg);
	}

	GMRFLib_openmp_implement_strategy(GMRFLib_OPENMP_PLACES_DEFAULT, NULL, NULL);

	return GMRFLib_SUCCESS;
}

GMRFLib_hgmrfm_type_tp GMRFLib_hgmrfm_what_type(int node, GMRFLib_hgmrfm_arg_tp *a)
{
	int i;
	GMRFLib_hgmrfm_type_tp t = { GMRFLib_HGMRFM_TP___VOID, -1, -1 };
	if (node < a->n + a->n_ext) {
		t.tp = GMRFLib_HGMRFM_TP_ETA;
		t.idx = node;
		t.tp_idx = 0;
	} else if (a->nf && node < a->idx_map_f[a->nf]) {
		t.tp = GMRFLib_HGMRFM_TP_F;
		for (i = 0; i < a->nf; i++) {
			if (node < a->idx_map_f[i + 1]) {
				t.tp_idx = i;
				t.idx = node - a->idx_map_f[i];
				break;
			}
		}
	} else if (a->nbeta && node < a->idx_map_beta[a->nbeta]) {
		t.tp = GMRFLib_HGMRFM_TP_BETA;
		t.tp_idx = node - a->idx_map_beta[0];
		t.idx = 0;
	} else {
		GMRFLib_ASSERT_RETVAL(0 == 1, GMRFLib_ESNH, t);
	}
	return t;
}

double GMRFLib_hgmrfm_Qfunc(int thread_id, int node, int nnode, double *UNUSED(values), void *arg)
{
	if (nnode < 0) {
		return NAN;
	}

	/*
	 * this is Qfunction for the hgmrfm-function 
	 */
	GMRFLib_hgmrfm_arg_tp *a = NULL;
	int ii, jj, equal;
	double value = 0.0;
	GMRFLib_hgmrfm_type_tp it, jt;

	a = (GMRFLib_hgmrfm_arg_tp *) arg;

	if (node < nnode) {
		ii = node;
		jj = nnode;
		equal = 0;
	} else {
		ii = nnode;
		jj = node;
		equal = (ii == jj);
	}

	// old non-caching code: it = GMRFLib_hgmrfm_what_type(ii, a); jt = GMRFLib_hgmrfm_what_type(jj, a);
	it = a->what_type[ii];
	jt = a->what_type[jj];

	if (equal || GMRFLib_graph_is_nb(ii, jj, a->eta_graph)) {
		value += a->eta_Q->Qfunc(thread_id, ii, jj, NULL, a->eta_Q->Qfunc_arg);
	}
	if (a->lc_Q && (equal || GMRFLib_graph_is_nb(ii, jj, a->lc_graph))) {
		value += a->lc_Q->Qfunc(thread_id, ii, jj, NULL, a->lc_Q->Qfunc_arg);
	}
	switch (it.tp) {
	case GMRFLib_HGMRFM_TP_ETA:
	{
		switch (jt.tp) {
		case GMRFLib_HGMRFM_TP_ETA:
			if (a->eta_ext_graph) {
				if (equal || GMRFLib_graph_is_nb(ii, jj, a->eta_ext_graph)) {
					value += a->eta_ext_Q->Qfunc(thread_id, ii, jj, NULL, a->eta_ext_Q->Qfunc_arg);
				}
			}
			return value;
		default:
			return value;
		}
		GMRFLib_ASSERT_RETVAL(0 == 1, GMRFLib_ESNH, 0.0);
	}
		break;

	case GMRFLib_HGMRFM_TP_F:
	{
		switch (jt.tp) {
		case GMRFLib_HGMRFM_TP_F:
			if (it.tp_idx == jt.tp_idx) {
				if ((it.idx == jt.idx) || GMRFLib_graph_is_nb(it.idx, jt.idx, a->f_graph[it.tp_idx])) {
					value += a->f_Qfunc[it.tp_idx] (thread_id, it.idx, jt.idx, NULL, a->f_Qfunc_arg[it.tp_idx]);
				}
			}
			/*
			 * only for the same index and different types; used to define `interaction between fields'. this is a 'workaround' for a INLA problem.. 
			 */
			if (a->ff_Qfunc) {
				if ((it.idx == jt.idx) && (it.tp_idx != jt.tp_idx) && a->ff_Qfunc[it.tp_idx][jt.tp_idx]) {
					value +=
					    a->ff_Qfunc[it.tp_idx][jt.tp_idx] (thread_id, it.idx, jt.idx, NULL,
									       a->ff_Qfunc_arg[it.tp_idx][jt.tp_idx]);
				}
			}
			return value;
		case GMRFLib_HGMRFM_TP_BETA:
			return value;
		case GMRFLib_HGMRFM_TP_LC:
			return value;
		default:
			GMRFLib_ASSERT_RETVAL(0 == 1, GMRFLib_ESNH, 0.0);
		}
		GMRFLib_ASSERT_RETVAL(0 == 1, GMRFLib_ESNH, 0.0);
	}
		break;

	case GMRFLib_HGMRFM_TP_BETA:
	{
		switch (jt.tp) {
		case GMRFLib_HGMRFM_TP_BETA:
			if (it.tp_idx == jt.tp_idx) {
				value += (a->prior_precision ? a->prior_precision[it.tp_idx] : 0.0);
			}
			return value;
		case GMRFLib_HGMRFM_TP_LC:
			return value;
		default:
			GMRFLib_ASSERT_RETVAL(0 == 1, GMRFLib_ESNH, 0.0);
		}
		GMRFLib_ASSERT_RETVAL(0 == 1, GMRFLib_ESNH, 0.0);
	}
		break;

	case GMRFLib_HGMRFM_TP_LC:
	{
		return value;
	}
		break;

	default:
		GMRFLib_ASSERT_RETVAL(0 == 1, GMRFLib_ESNH, 0.0);
		break;
	}

	return value;
}

int GMRFLib_free_hgmrfm(GMRFLib_hgmrfm_tp *hgmrfm)
{
	if (!hgmrfm) {
		return GMRFLib_SUCCESS;
	}
	GMRFLib_hgmrfm_arg_tp *a = (GMRFLib_hgmrfm_arg_tp *) hgmrfm->Qfunc_arg;

	GMRFLib_graph_free(hgmrfm->graph);
	GMRFLib_free_constr(hgmrfm->constr);
	GMRFLib_free_tabulate_Qfunc(a->eta_Q);
	GMRFLib_graph_free(a->eta_graph);
	GMRFLib_free_tabulate_Qfunc(a->lc_Q);
	GMRFLib_graph_free(a->lc_graph);
	Free(a->idx_map_f);
	Free(a->idx_map_beta);
	Free(a->idx_map_lc);
	Free(a->what_type);
	Free(a);
	Free(hgmrfm);

	return GMRFLib_SUCCESS;
}

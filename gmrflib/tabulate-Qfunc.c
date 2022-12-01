
/* tabulate-Qfunc.c
 * 
 * Copyright (C) 2004-2022 Havard Rue
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

#include <stdio.h>
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif
#include <stdlib.h>

#ifndef GITCOMMIT
#define GITCOMMIT
#endif
static const char GitID[] = "file: " __FILE__ "  " GITCOMMIT;

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

#define TAB_FUNC_CORE(_prec_scale)					\
	GMRFLib_tabulate_Qfunc_arg_tp *args = NULL;			\
	double prec = 1.0;						\
	args = (GMRFLib_tabulate_Qfunc_arg_tp *) arg;			\
	if (_prec_scale) {						\
		prec = GMRFLib_SET_PREC(args);				\
	}								\
	if (nnode >= 0) {						\
		int imin, imax;						\
		if (node <= nnode) {					\
			imin = node;					\
			imax = nnode;					\
		} else {						\
			imin = nnode;					\
			imax = node;					\
		}							\
		double *dp = NULL;					\
		if (args->Q) {						\
			int offset = args->Q->s->ia[imin];		\
			int j = offset + GMRFLib_iwhich_sorted(imax, offset + args->Q->s->ja, args->Q->s->ia[imin + 1] - offset); \
			dp = &(args->Q->a[j]);				\
		} else if (args->Q_idx) {				\
			int ii = -1;					\
			map_ii_get(args->Q_idx[imin], imax, &ii);	\
			dp = &(args->Q->a[ii]);				\
		} else {						\
			dp = map_id_ptr(args->values[imin], imax);	\
		}							\
									\
		if (_prec_scale) {					\
			val = prec * (*dp);				\
		} else {						\
			val = *dp;					\
		}							\
	} else {							\
		int len = 0;						\
		if (args->Q) {						\
			int j = args->Q->s->ia[node];			\
			len = args->Q->s->ia[node + 1] - j;		\
			Memcpy(values, &(args->Q->a[j]), len * sizeof(double)); \
		} else {						\
			val = NAN;					\
		}							\
									\
		if (_prec_scale) {					\
			GMRFLib_dscale(len, prec, values);			\
		}							\
	}

double GMRFLib_tabulate_Qfunction(int thread_id, int node, int nnode, double *values, void *arg)
{
	double val = 0.0;
	TAB_FUNC_CORE(1);
	return val;
}

double GMRFLib_tabulate_Qfunction_std(int thread_id, int node, int nnode, double *values, void *arg)
{
	double val = 0.0;
	TAB_FUNC_CORE(0);
	return val;
}

int GMRFLib_tabulate_Qfunc(int thread_id,
			   GMRFLib_tabulate_Qfunc_tp ** tabulate_Qfunc, GMRFLib_graph_tp * graph,
			   GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_arg, double **log_prec_omp)
{
	return (GMRFLib_tabulate_Qfunc_core(thread_id, tabulate_Qfunc, graph, Qfunc, Qfunc_arg, log_prec_omp, 0));
}

int GMRFLib_tabulate_Qfunc_core(int thread_id,
				GMRFLib_tabulate_Qfunc_tp ** tabulate_Qfunc, GMRFLib_graph_tp * graph,
				GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_arg, double **log_prec_omp, int force)
{
	int i, j, k;
	*tabulate_Qfunc = Calloc(1, GMRFLib_tabulate_Qfunc_tp);

	if (!force) {
		// 
		// this is VERY SPECIAL: with this we can reuse code with the preopt functionality
		// 
		if (Qfunc == GMRFLib_preopt_Qfunc) {
			(*tabulate_Qfunc)->Qfunc = Qfunc;
			(*tabulate_Qfunc)->Qfunc_arg = Qfunc_arg;
			return GMRFLib_SUCCESS;
		}
	}

	GMRFLib_tabulate_Qfunc_arg_tp *arg = NULL;
	arg = Calloc(1, GMRFLib_tabulate_Qfunc_arg_tp);
	(*tabulate_Qfunc)->Qfunc_arg = (void *) arg;

	arg->n = graph->n;
	if (log_prec_omp) {
		int tmax = GMRFLib_MAX_THREADS();
		arg->log_prec_omp = Calloc(tmax, double *);
		for (i = 0; i < tmax; i++) {
			arg->log_prec_omp[i] = log_prec_omp[i];
		}
		(*tabulate_Qfunc)->Qfunc = GMRFLib_tabulate_Qfunction;
	} else {
		arg->log_prec_omp = NULL;
		(*tabulate_Qfunc)->Qfunc = GMRFLib_tabulate_Qfunction_std;
	}

	if (GMRFLib_smtp == GMRFLib_SMTP_PARDISO || GMRFLib_smtp == GMRFLib_SMTP_TAUCS) {
		GMRFLib_Q2csr(thread_id, &(arg->Q), graph, Qfunc, Qfunc_arg);
		if (arg->Q->a[0] < 0.0 || ISNAN(arg->Q->a[0]) || ISINF(arg->Q->a[0])) {
			P(arg->Q->a[0]);
		}
		assert(arg->Q->a[0] >= 0.0);
		GMRFLib_graph_duplicate(&(arg->graph), graph);
	} else {
		arg->values = Calloc(graph->n, map_id *);
		map_id *work = Calloc(graph->n, map_id);
		for (i = 0; i < graph->n; i++) {
			arg->values[i] = work + i;
		}

		omp_set_num_threads(GMRFLib_openmp->max_threads_inner);
//#pragma omp parallel for private(i, j, k) num_threads(GMRFLib_openmp->max_threads_inner)
		for (i = 0; i < graph->n; i++) {
			map_id_init_hint(arg->values[i], graph->lnnbs[i] + 1);
			map_id_set(arg->values[i], i, (*Qfunc) (thread_id, i, i, NULL, Qfunc_arg));	/* diagonal */

			for (j = 0; j < graph->lnnbs[i]; j++) {
				k = graph->lnbs[i][j];
				map_id_set(arg->values[i], k, (*Qfunc) (thread_id, i, k, NULL, Qfunc_arg));
			}
		}
	}

	return GMRFLib_SUCCESS;
}

int GMRFLib_tabulate_Qfunc_from_file(GMRFLib_tabulate_Qfunc_tp ** tabulate_Qfunc, GMRFLib_graph_tp ** graph, const char *filename,
				     int dim, double **log_prec_omp)
{
	/*
	 * as GMRFLib_tabulate_Qfunc(), but reads the Q_ij values from file with name FILENAME, in format
	 * 
	 * i j Q_{ij} : : : i j Q_{ij}
	 *
	 *
	 * only values i<=j are required. duplicated values are added up.
	 */

	const int debug = 0;
	int i, j, ii, jj, k, ntriples, err, imin = INT_MAX, jmin = INT_MAX, off = 0, sparse = 0;
	double value;

	GMRFLib_tabulate_Qfunc_arg_tp *arg = NULL;
	GMRFLib_io_tp *io = NULL;
	GMRFLib_error_handler_tp *old_handler;
	GMRFLib_matrix_tp *M = NULL;

	/*
	 * step 1. build the graph 
	 */

	GMRFLib_ged_tp *ged = NULL;

	GMRFLib_ged_init(&ged, NULL);

	/*
	 * to fix the dimension, possibly padding with zero's 
	 */
	GMRFLib_ged_add(ged, 0, 0);
	if (dim > 0) {
		GMRFLib_ged_add(ged, dim - 1, dim - 1);
	}

	if (GMRFLib_is_fmesher_file(filename, (long int) 0, -1) == GMRFLib_SUCCESS) {
		M = GMRFLib_read_fmesher_file(filename, (long int) 0, -1);
		sparse = (M->i && M->j);
		if (!sparse) {
			assert(M->ncol == 3);
		} else {
			/*
			 * make sure to fix the dimension of the matrix
			 */
			assert(M->nrow == M->ncol);
			GMRFLib_ged_add(ged, M->nrow - 1, M->ncol - 1);
		}
		if (sparse) {
			for (k = 0; k < M->elems; k++) {
				i = M->i[k];
				j = M->j[k];
				imin = IMIN(imin, i);
				jmin = IMIN(jmin, j);
			}
		} else {
			for (k = 0; k < M->nrow; k++) {
				i = (int) M->A[k + 0 * M->nrow];
				j = (int) M->A[k + 1 * M->nrow];
				imin = IMIN(imin, i);
				jmin = IMIN(jmin, j);
			}
		}

		if (sparse) {
			ntriples = M->elems;
			for (k = 0; k < M->elems; k++) {
				i = M->i[k];
				j = M->j[k];
				value = M->values[k];
				if (debug) {
					printf("read (i,j,val) = (%d,%d,%g)\n", i, j, value);
				}
				GMRFLib_ged_add(ged, i, j);
			}
		} else {
			ntriples = M->nrow;
			for (k = 0; k < M->nrow; k++) {
				i = (int) M->A[k + 0 * M->nrow];
				j = (int) M->A[k + 1 * M->nrow];
				value = M->A[k + 2 * M->nrow];
				if (debug) {
					printf("read (i,j,val) = (%d,%d,%g)\n", i, j, value);
				}
				GMRFLib_ged_add(ged, i, j);
			}
		}
		/*
		 * I will free matrix M later... 
		 */
	} else {
		/*
		 * read it first to determine if this is a zero-based or one-based graph 
		 */
		GMRFLib_EWRAP0(GMRFLib_io_open(&io, filename, "r"));
		while (1) {
			old_handler = GMRFLib_set_error_handler_off();
			err = GMRFLib_io_read_next(io, &i, "%d");
			GMRFLib_set_error_handler(old_handler);

			if (err == GMRFLib_SUCCESS) {
				/*
				 * then the rest must be present to 
				 */
				GMRFLib_EWRAP0(GMRFLib_io_read_next(io, &j, "%d"));
				GMRFLib_EWRAP0(GMRFLib_io_read_next(io, &value, "%lf"));
				imin = IMIN(imin, i);
				jmin = IMIN(jmin, j);
			} else {
				break;
			}
		}
		GMRFLib_EWRAP0(GMRFLib_io_close(io));

		GMRFLib_ASSERT(((imin == 0 || imin == 1) && (jmin == 0 || jmin == 1)), GMRFLib_ESNH);
		off = (IMIN(imin, jmin) == 1 ? 1 : 0);

		ntriples = 0;
		GMRFLib_EWRAP0(GMRFLib_io_open(&io, filename, "r"));
		while (1) {
			old_handler = GMRFLib_set_error_handler_off();
			err = GMRFLib_io_read_next(io, &i, "%d");
			GMRFLib_set_error_handler(old_handler);

			if (err == GMRFLib_SUCCESS) {
				/*
				 * then the rest must be present to 
				 */
				GMRFLib_EWRAP0(GMRFLib_io_read_next(io, &j, "%d"));
				GMRFLib_EWRAP0(GMRFLib_io_read_next(io, &value, "%lf"));
				GMRFLib_ged_add(ged, i - off, j - off);

				if (debug)
					printf("read (i,j,val) = (%d,%d,%g)\n", i, j, value);
				ntriples++;
			} else {
				break;
			}
		}
		GMRFLib_EWRAP0(GMRFLib_io_close(io));
	}

	GMRFLib_ged_build(graph, ged);
	GMRFLib_ged_free(ged);

	/*
	 * step 2. build the tabulate_Qfunc structure 
	 */
	*tabulate_Qfunc = Calloc(1, GMRFLib_tabulate_Qfunc_tp);
	(*tabulate_Qfunc)->Qfunc = GMRFLib_tabulate_Qfunction; /* the Qfunction to use */
	arg = Calloc(1, GMRFLib_tabulate_Qfunc_arg_tp);
	(*tabulate_Qfunc)->Qfunc_arg = (void *) arg;

	arg->n = (*graph)->n;
	arg->values = Calloc((*graph)->n, map_id *);
	if (log_prec_omp) {
		int tmax = GMRFLib_MAX_THREADS();
		arg->log_prec_omp = Calloc(tmax, double *);
		for (i = 0; i < tmax; i++) {
			arg->log_prec_omp[i] = log_prec_omp[i];
		}
	} else {
		arg->log_prec_omp = NULL;
	}

	/*
	 * allocate hash-table with the *correct* number of elements
	 */
	map_id *work = Calloc((*graph)->n, map_id);
	for (i = 0; i < (*graph)->n; i++) {
		arg->values[i] = work + i;
	}
	for (i = 0; i < (*graph)->n; i++) {
		map_id_init_hint(arg->values[i], (*graph)->lnnbs[i] + 1);
	}

	/*
	 * fill all entries in the graph with zero 
	 */
	for (i = 0; i < (*graph)->n; i++) {
		map_id_set(arg->values[i], i, 0.0);
		for (jj = 0; jj < (*graph)->lnnbs[i]; jj++) {
			j = (*graph)->lnbs[i][jj];
			map_id_set(arg->values[i], j, 0.0);
		}
	}

	/*
	 * then read the file again 
	 */
	if (M) {
		/*
		 * ...or we have it already 
		 */
		if (sparse) {
			for (k = 0; k < M->elems; k++) {
				i = M->i[k];
				j = M->j[k];
				value = M->values[k];
				if (i <= j) {
					i = i - off;
					j = j - off;
					ii = IMIN(i, j);
					jj = IMAX(i, j);
					map_id_set(arg->values[ii], jj, value);
					if (debug) {
						printf("set (i,j,val) = (%d,%d,%g)\n", i, j, value);
					}
				}
			}
		} else {
			for (k = 0; k < M->nrow; k++) {
				i = (int) M->A[k + 0 * M->nrow];
				j = (int) M->A[k + 1 * M->nrow];
				value = M->A[k + 2 * M->nrow];

				if (i <= j) {
					i = i - off;
					j = j - off;
					ii = IMIN(i, j);
					jj = IMAX(i, j);
					map_id_set(arg->values[ii], jj, value);
					if (debug) {
						printf("set (i,j,val) = (%d,%d,%g)\n", i, j, value);
					}
				}
			}
		}
	} else {
		GMRFLib_EWRAP0(GMRFLib_io_open(&io, filename, "r"));
		k = 0;
		while (k < ntriples) {
			GMRFLib_EWRAP0(GMRFLib_io_read_next(io, &i, "%d"));
			GMRFLib_EWRAP0(GMRFLib_io_read_next(io, &j, "%d"));
			GMRFLib_EWRAP0(GMRFLib_io_read_next(io, &value, "%lf"));

			if (i <= j) {
				i = i - off;
				j = j - off;
				ii = IMIN(i, j);
				jj = IMAX(i, j);
				map_id_set(arg->values[ii], jj, value);
				if (debug)
					printf("set (i,j,val) = (%d,%d,%g)\n", i, j, value);
			}
			k++;
		}
		GMRFLib_EWRAP0(GMRFLib_io_close(io));
	}

	GMRFLib_matrix_free(M);
	return GMRFLib_SUCCESS;
}

int GMRFLib_tabulate_Qfunc_from_list(GMRFLib_tabulate_Qfunc_tp ** tabulate_Qfunc, GMRFLib_graph_tp ** graph,
				     int ntriples, int *ilist, int *jlist, double *Qijlist, int dim, double **log_prec_omp)
{
	/*
	 * as GMRFLib_tabulate_Qfunc(), but get its Q_ij values from its arguments
	 * 
	 * i j Q_{ij} : : : i j Q_{ij}
	 * 
	 */

	int i, imin = INT_MAX, jmin = INT_MAX, off;
	GMRFLib_tabulate_Qfunc_arg_tp *arg = NULL;

	/*
	 * step 1. build the graph 
	 */
	GMRFLib_ged_tp *ged = NULL;

	GMRFLib_ged_init(&ged, NULL);

	for (i = 0; i < ntriples; i++) {
		imin = IMIN(imin, ilist[i]);
		jmin = IMIN(jmin, jlist[i]);
	}
	GMRFLib_ASSERT(((imin == 0 || imin == 1) && (jmin == 0 || jmin == 1)), GMRFLib_ESNH);
	off = (IMIN(imin, jmin) == 1 ? 1 : 0);

	/*
	 * to fix the dimension, possibly padding with zero's 
	 */
	GMRFLib_ged_add(ged, 0, 0);
	if (dim > 0) {
		GMRFLib_ged_add(ged, dim - 1, dim - 1);
	}

	for (i = 0; i < ntriples; i++) {
		GMRFLib_ged_add(ged, ilist[i] - off, jlist[i] - off);
	}

	GMRFLib_ged_build(graph, ged);
	GMRFLib_graph_prepare(*graph);
	GMRFLib_ged_free(ged);

	/*
	 * step 2. build the tabulate_Qfunc structure 
	 */
	*tabulate_Qfunc = Calloc(1, GMRFLib_tabulate_Qfunc_tp);
	(*tabulate_Qfunc)->Qfunc = GMRFLib_tabulate_Qfunction; /* the Qfunction to use */
	arg = Calloc(1, GMRFLib_tabulate_Qfunc_arg_tp);
	(*tabulate_Qfunc)->Qfunc_arg = (void *) arg;

	arg->n = (*graph)->n;
	arg->values = Calloc((*graph)->n, map_id *);
	if (log_prec_omp) {
		int tmax = GMRFLib_MAX_THREADS();
		arg->log_prec_omp = Calloc(tmax, double *);
		for (i = 0; i < tmax; i++) {
			arg->log_prec_omp[i] = log_prec_omp[i];
		}
	} else {
		arg->log_prec_omp = NULL;
	}

	/*
	 * allocate hash-table with the *correct* number of elements
	 */
	map_id *work = Calloc((*graph)->n, map_id);
	for (i = 0; i < (*graph)->n; i++) {
		arg->values[i] = work + i;
	}

//#pragma omp parallel for private(i)
	for (i = 0; i < (*graph)->n; i++) {
		int j, jj;
		map_id_init_hint(arg->values[i], (*graph)->lnnbs[i] + 1);
		map_id_set(arg->values[i], i, 0.0);
		for (jj = 0; jj < (*graph)->lnnbs[i]; jj++) {
			j = (*graph)->lnbs[i][jj];
			map_id_set(arg->values[i], j, 0.0);    /* fill them with default = 0.0 */
		}
	}

	for (i = 0; i < ntriples; i++) {
		int ii, jj;
		if (ilist[i] <= jlist[i]) {
			ii = ilist[i] - off;
			jj = jlist[i] - off;
			map_id_set(arg->values[ii], jj, Qijlist[i]);
		}
	}

	return GMRFLib_SUCCESS;
}

int GMRFLib_tabulate_Qfunc_from_list2(GMRFLib_tabulate_Qfunc_tp ** tabulate_Qfunc, GMRFLib_graph_tp * graph,
				      int ntriples, int *ilist, int *jlist, double *Qijlist, int UNUSED(dim), double **log_prec_omp)
{
	// this is a special version for Qfunc_rgeneric, as we assume here that graph is know.

	/*
	 * as GMRFLib_tabulate_Qfunc(), but get its Q_ij values from its arguments
	 * 
	 * i j Q_{ij} : : : i j Q_{ij}
	 * 
	 */

	int i, imin = INT_MAX, jmin = INT_MAX, off;
	GMRFLib_tabulate_Qfunc_arg_tp *arg = NULL;

	for (i = 0; i < ntriples; i++) {
		imin = IMIN(imin, ilist[i]);
		jmin = IMIN(jmin, jlist[i]);
	}
	GMRFLib_ASSERT(((imin == 0 || imin == 1) && (jmin == 0 || jmin == 1)), GMRFLib_ESNH);
	off = (IMIN(imin, jmin) == 1 ? 1 : 0);

	*tabulate_Qfunc = Calloc(1, GMRFLib_tabulate_Qfunc_tp);
	(*tabulate_Qfunc)->Qfunc = GMRFLib_tabulate_Qfunction; /* the Qfunction to use */
	arg = Calloc(1, GMRFLib_tabulate_Qfunc_arg_tp);
	(*tabulate_Qfunc)->Qfunc_arg = (void *) arg;

	arg->n = graph->n;
	arg->values = Calloc(graph->n, map_id *);
	if (log_prec_omp != NULL) {
		int tmax = GMRFLib_MAX_THREADS();
		arg->log_prec_omp = Calloc(tmax, double *);
		for (i = 0; i < tmax; i++) {
			arg->log_prec_omp[i] = log_prec_omp[i];
		}
	} else {
		arg->log_prec_omp = NULL;
	}

	/*
	 * allocate hash-table with the *correct* number of elements
	 */
	map_id *work = Calloc(graph->n, map_id);
	for (i = 0; i < graph->n; i++) {
		arg->values[i] = work + i;
	}
	for (i = 0; i < graph->n; i++) {
		int j, jj;

		map_id_init_hint(arg->values[i], graph->lnnbs[i] + 1);
		map_id_set(arg->values[i], i, 0.0);
		for (jj = 0; jj < graph->lnnbs[i]; jj++) {
			j = graph->lnbs[i][jj];
			map_id_set(arg->values[i], j, 0.0);    /* fill them with default = 0.0 */
		}
	}

	for (i = 0; i < ntriples; i++) {
		int ii, jj;
		if (ilist[i] <= jlist[i]) {
			ii = ilist[i] - off;
			jj = jlist[i] - off;
			map_id_set(arg->values[ii], jj, Qijlist[i]);
		}
	}

	return GMRFLib_SUCCESS;
}

int GMRFLib_free_tabulate_Qfunc(GMRFLib_tabulate_Qfunc_tp * tabulate_Qfunc)
{
	if (!tabulate_Qfunc)
		return GMRFLib_SUCCESS;
	// 
	// this is VERY SPECIAL: with this we can reuse code with the preopt functionality
	// 
	if (tabulate_Qfunc->Qfunc == GMRFLib_preopt_Qfunc) {
		Free(tabulate_Qfunc);
		return GMRFLib_SUCCESS;
	}

	int i;
	GMRFLib_tabulate_Qfunc_arg_tp *arg = NULL;

	if (tabulate_Qfunc) {
		arg = (GMRFLib_tabulate_Qfunc_arg_tp *) tabulate_Qfunc->Qfunc_arg;
		if (arg->graph) {
			GMRFLib_graph_free(arg->graph);
		}
		if (arg->Q) {
			GMRFLib_csr_free(&(arg->Q));
		}
		if (arg->Q_idx) {
			for (i = 0; i < arg->n; i++) {
				map_ii_free(arg->Q_idx[i]);
				Free(arg->Q_idx[i]);
			}
			Free(arg->Q_idx);
		}
		if (arg->values) {
			for (i = 0; i < arg->n; i++) {
				map_id_free(arg->values[i]);
			}
			Free(arg->values[0]);
			Free(arg->values);
		}

		Free(arg->log_prec_omp);
		Free(arg);
		Free(tabulate_Qfunc);
	}

	return GMRFLib_SUCCESS;
}


/* GMRFLib-smtp-taucs.c
 * 
 * Copyright (C) 2001-2006 Havard Rue
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

/*!
  \file smtp-taucs.c
  \brief The interface towards the TAUCS-library
*/

#include <stddef.h>
#include <assert.h>
#include <strings.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"
#include "amd.h"
#include "metis.h"

#ifndef HGVERSION
#define HGVERSION
#endif
static const char RCSId[] = "file: " __FILE__ "  " HGVERSION;

/* Pre-hg-Id: $Id: smtp-taucs.c,v 1.162 2010/02/27 08:32:38 hrue Exp $ */

/* 
   if TRUE, then we use my modified routine to convert from supernodal_factor to ccs, which preserves zeros in L. this gives
   speedup for the computations in Qinv. So far I have tested this works fine and is correct. So there two instances in the code
   where we make use of this; see below. 
 */
static int include_zeros_in_L = GMRFLib_TRUE;

/* 
   how large should `nset' be before doing `memset' instead of a `for' loop.
*/
#define GMRFLib_NSET_LIMIT(nset, size, n)  IMAX(10, (n)/10/(size))

/* 
   First some modified code from the TAUCS library. Checked to be ok version 2.0 and 2.2

   taucs_datatype is set to double in GMRFLib/taucs.h
*/
taucs_ccs_matrix *my_taucs_dsupernodal_factor_to_ccs(void *vL)
{
	/*
	 * this is to be called for a lower triangular double matrix only.
	 * 
	 * it includes also zero terms as long as i>=j.
	 * 
	 */

	supernodal_factor_matrix *L = (supernodal_factor_matrix *) vL;
	taucs_ccs_matrix *C = NULL;
	int n, nnz;
	int i, j, ip, jp, sn, next;
	taucs_datatype v;
	int *len = NULL;

	n = L->n;

	len = Malloc(n, int);

	if (!len) {
		return NULL;
	}
	nnz = 0;
	for (sn = 0; sn < L->n_sn; sn++) {
		for (jp = 0; jp < L->sn_size[sn]; jp++) {
			j = L->sn_struct[sn][jp];
			len[j] = 0;

			for (ip = jp; ip < L->sn_size[sn]; ip++) {
				i = L->sn_struct[sn][ip];
				if (i >= j) {
					len[j]++;
					nnz++;
				}
			}
			for (ip = L->sn_size[sn]; ip < L->sn_up_size[sn]; ip++) {
				i = L->sn_struct[sn][ip];
				if (i >= j) {
					len[j]++;
					nnz++;
				}
			}
		}
	}

	C = taucs_dccs_create(n, n, nnz);
	if (!C) {
		free(len);
		return NULL;
	}
	C->flags = TAUCS_DOUBLE;
	C->flags |= TAUCS_TRIANGULAR | TAUCS_LOWER;	       /* this was a bug in version 2.0 of taucs */

	(C->colptr)[0] = 0;
	for (j = 1; j <= n; j++) {
		(C->colptr)[j] = (C->colptr)[j - 1] + len[j - 1];
	}

	free(len);
	for (sn = 0; sn < L->n_sn; sn++) {
		for (jp = 0; jp < L->sn_size[sn]; jp++) {
			j = L->sn_struct[sn][jp];
			next = (C->colptr)[j];

			for (ip = jp; ip < L->sn_size[sn]; ip++) {
				i = L->sn_struct[sn][ip];
				v = L->sn_blocks[sn][jp * L->sn_blocks_ld[sn] + ip];

				if (i >= j) {
					(C->rowind)[next] = i;
					(C->values.d)[next] = v;
					next++;
				}
			}
			for (ip = L->sn_size[sn]; ip < L->sn_up_size[sn]; ip++) {
				i = L->sn_struct[sn][ip];
				v = L->up_blocks[sn][jp * L->up_blocks_ld[sn] + (ip - L->sn_size[sn])];

				if (i >= j) {
					(C->rowind)[next] = i;
					(C->values.d)[next] = v;
					next++;
				}
			}
		}
	}
	return C;
}

/* 
   copy a supernodal_factor_matrix
*/
supernodal_factor_matrix *GMRFLib_sm_fact_duplicate_TAUCS(supernodal_factor_matrix * L)
{
#define DUPLICATE(name,len,type) if (1) {					\
		if (L->name && ((len) > 0)) {				\
			LL->name = (type *)Calloc((len), type);		\
			memcpy(LL->name,L->name,(size_t)(len)*sizeof(type)); \
		} else {						\
			LL->name = (type *)NULL;			\
		}							\
	}

	supernodal_factor_matrix *LL = NULL;
	int n, np, n_sn;

	if (!L) {
		return NULL;
	}
	LL = (supernodal_factor_matrix *) Calloc((size_t) 1, supernodal_factor_matrix);
	LL->flags = L->flags;
	LL->uplo = L->uplo;
	LL->n = L->n;
	LL->n_sn = L->n_sn;

	n = LL->n;
	np = n + 1;
	n_sn = LL->n_sn;

	DUPLICATE(sn_size, np, int);
	DUPLICATE(sn_up_size, np, int);
	DUPLICATE(first_child, np, int);
	DUPLICATE(next_child, np, int);
	DUPLICATE(parent, np, int);
	DUPLICATE(sn_blocks_ld, n_sn, int);
	DUPLICATE(up_blocks_ld, n_sn, int);

	{
		{
			int i;
			LL->sn_struct = (int **) Calloc(n_sn, int *);
			for (i = 0; i < LL->n_sn; i++) {
				DUPLICATE(sn_struct[i], LL->sn_up_size[i], int);
			}
		}
		{
			int i;
			LL->sn_blocks = (double **) Calloc(n_sn, double *);
			for (i = 0; i < LL->n_sn; i++) {
				DUPLICATE(sn_blocks[i], ISQR(LL->sn_size[i]), double);
			}
		}
		{
			int i;
			LL->up_blocks = (double **) Calloc(n_sn, double *);
			for (i = 0; i < LL->n_sn; i++) {
				DUPLICATE(up_blocks[i], (LL->sn_up_size[i] - LL->sn_size[i]) * (LL->sn_size)[i], double);
			}
		}
	}

#undef DUPLICATE
	return LL;
}

void taucs_ccs_metis5(taucs_ccs_matrix * m, int **perm, int **invperm, char *which)
{
	// this for metis version 5

	int n, nnz, i, j, ip;
	int *xadj;
	int *adj;
	int *len;
	int *ptr;
	int ret;

	assert(sizeof(idx_t) == sizeof(int));

	if (!(m->flags & TAUCS_SYMMETRIC) && !(m->flags & TAUCS_HERMITIAN)) {
		taucs_printf("taucs_ccs_treeorder: METIS ordering only works on symmetric matrices.\n");
		*perm = NULL;
		*invperm = NULL;
		return;
	}
	/*
	 * this routine may actually work on UPPER as well 
	 */
	if (!(m->flags & TAUCS_LOWER)) {
		taucs_printf("taucs_ccs_metis: the lower part of the matrix must be represented.\n");
		*perm = NULL;
		*invperm = NULL;
		return;
	}

	n = m->n;
	nnz = (m->colptr)[n];

	*perm = Calloc(n, int);
	*invperm = Calloc(n, int);

	xadj = Calloc(n + 1, int);
	adj = Calloc(2 * nnz, int);


	if (!(*perm) || !(*invperm) || !xadj || !adj) {
		Free(*perm);
		Free(*invperm);
		Free(xadj);
		Free(adj);
		*perm = *invperm = NULL;
		return;
	}

	ptr = len = *perm;
	for (i = 0; i < n; i++)
		len[i] = 0;

	for (j = 0; j < n; j++) {
		for (ip = (m->colptr)[j]; ip < (m->colptr)[j + 1]; ip++) {
			/*
			 * i = (m->rowind)[ip] - (m->indshift);
			 */
			i = (m->rowind)[ip];
			if (i != j) {
				len[i]++;
				len[j]++;
			}
		}
	}

	xadj[0] = 0;
	for (i = 1; i <= n; i++)
		xadj[i] = xadj[i - 1] + len[i - 1];

	for (i = 0; i < n; i++)
		ptr[i] = xadj[i];

	for (j = 0; j < n; j++) {
		for (ip = (m->colptr)[j]; ip < (m->colptr)[j + 1]; ip++) {
			/*
			 * i = (m->rowind)[ip] - (m->indshift);
			 */
			i = (m->rowind)[ip];
			if (i != j) {
				adj[ptr[i]] = j;
				adj[ptr[j]] = i;
				ptr[i]++;
				ptr[j]++;
			}
		}
	}
	idx_t options[METIS_NOPTIONS];
	// Have to adapt to the PARDISO metis libs
	// METIS_SetDefaultOptions(options);
	for (i = 0; i < METIS_NOPTIONS; i++) {
		options[i] = -1;
	}

	options[METIS_OPTION_NUMBERING] = 0;
	options[METIS_OPTION_NSEPS] = 5;
	options[METIS_OPTION_COMPRESS] = 0;
	options[METIS_OPTION_PFACTOR] = 100;

#if defined(NO_PARDISO_LIB)
	// this the metis5 lib
	ret = METIS_NodeND(&n, xadj, adj, NULL, options, *perm, *invperm);
#else
	// this is defined in the pardiso libs
	ret = METIS51_NodeND(&n, xadj, adj, NULL, options, *perm, *invperm);
#endif
	if (ret != METIS_OK)
		GMRFLib_ERROR(GMRFLib_EREORDER);

	Free(xadj);
	Free(adj);
}

GMRFLib_sizeof_tp GMRFLib_sm_fact_nnz_TAUCS(supernodal_factor_matrix * L)
{
	/*
	 * return the number of non-zeros in the matrix 
	 */
	GMRFLib_sizeof_tp nnz = 0;
	int jp, sn;

	for (sn = 0; sn < L->n_sn; sn++) {
		for (jp = 0; jp < L->sn_size[sn]; jp++) {
			nnz += L->sn_size[sn] - jp;
			nnz += L->sn_up_size[sn] - L->sn_size[sn];
		}
	}
	return (nnz);
}

GMRFLib_sizeof_tp GMRFLib_sm_fact_sizeof_TAUCS(supernodal_factor_matrix * L)
{
	/*
	 * return, approximately, the size of L 
	 */
	GMRFLib_sizeof_tp siz = 0;

	if (!L) {
		return siz;
	}

	int n, np, n_sn, i;

	siz += sizeof(supernodal_factor_matrix);
	n = L->n;
	np = n + 1;
	n_sn = L->n_sn;
	siz += np * sizeof(int) * 5 + n_sn * sizeof(int) * 2;
	siz += n * sizeof(int *);
	if (L->sn_up_size) {
		for (i = 0; i < L->n_sn; i++) {
			siz += L->sn_up_size[i] * sizeof(int);
		}
	}

	siz += n_sn * sizeof(double *);
	if (L->sn_size) {
		for (i = 0; i < L->n_sn; i++) {
			siz += ISQR(L->sn_size[i]) * sizeof(double);
		}
	}

	siz += n_sn * sizeof(double *);
	if (L->sn_up_size && L->sn_size) {
		for (i = 0; i < L->n_sn; i++) {
			siz += (L->sn_up_size[i] - L->sn_size[i]) * (L->sn_size)[i] * sizeof(double);
		}
	}
	return siz;
}

/* 
   make a copy of a ccs-matrix
*/
taucs_ccs_matrix *GMRFLib_L_duplicate_TAUCS(taucs_ccs_matrix * L, int flags)
{
	/*
	 * copy a square matrix 
	 */

	taucs_ccs_matrix *LL = NULL;
	int n, nnz;

	if (!L) {
		return NULL;
	}

	n = L->n;
	nnz = L->colptr[L->n];
	LL = taucs_ccs_create(n, n, nnz, flags);

	memcpy(LL->colptr, L->colptr, (n + 1) * sizeof(int));
	memcpy(LL->rowind, L->rowind, nnz * sizeof(int));
	memcpy(LL->values.d, L->values.d, nnz * sizeof(double));

	return LL;
}

int GMRFLib_print_ccs_matrix(FILE * fp, taucs_ccs_matrix * L)
{
	if (!L) {
		return GMRFLib_SUCCESS;
	}

	int i;
	int n = L->n;
	int nnz = L->colptr[L->n];

	fprintf(fp, "n = %d\n", n);
	fprintf(fp, "nnz = %d\n", nnz);

	for (i = 0; i < n + 1; i++) {
		fprintf(fp, "\tcolptr[%1d] = %1d\n", i, L->colptr[i]);
	}
	for (i = 0; i < nnz; i++) {
		fprintf(fp, "\trowind[%1d] = %1d\n", i, L->rowind[i]);
		fprintf(fp, "\tvalues[%1d] = %.12g\n", i, L->values.d[i]);
	}

	return GMRFLib_SUCCESS;
}

GMRFLib_sizeof_tp GMRFLib_L_sizeof_TAUCS(taucs_ccs_matrix * L)
{
	/*
	 * return, approximately, the sizeof L 
	 */

	GMRFLib_sizeof_tp siz = 0;
	int n, nnz;

	if (!L) {
		return siz;
	}

	n = L->n;
	nnz = L->colptr[L->n];
	siz += (n + 1) * sizeof(int) + nnz * sizeof(int) + nnz * sizeof(double);

	return siz;
}

int GMRFLib_compute_reordering_TAUCS(int **remap, GMRFLib_graph_tp * graph, GMRFLib_reorder_tp reorder, GMRFLib_global_node_tp * gn_ptr)
{
	/*
	 * new improved version which treats global nodes spesifically. 
	 */
	int i, j, k, ic, ne, n, ns, nnz, *perm = NULL, *iperm = NULL, limit, free_subgraph, *iperm_new = NULL, simple;
	char *fixed = NULL, *p = NULL;
	taucs_ccs_matrix *Q = NULL;
	GMRFLib_graph_tp *subgraph = NULL;

	if (!graph || graph->n == 0) {
		return GMRFLib_SUCCESS;
	}

	if (reorder == GMRFLib_REORDER_IDENTITY || reorder == GMRFLib_REORDER_REVERSE_IDENTITY) {
		int *imap = Calloc(graph->n, int);
		if (reorder == GMRFLib_REORDER_IDENTITY) {
			for (i = 0; i < graph->n; i++) {
				imap[i] = i;
			}
		} else if (reorder == GMRFLib_REORDER_REVERSE_IDENTITY) {
			for (i = 0; i < graph->n; i++) {
				imap[i] = graph->n - 1 - i;
			}
		} else {
			assert(0 == 1);
		}
		*remap = imap;
		return GMRFLib_SUCCESS;
	}

	/*
	 * check if we have a simple solution --> no neigbours 
	 */
	for (i = 0, simple = 1; i < graph->n && simple; i++) {
		simple = (graph->nnbs[i] > 0 ? 0 : 1);
	}
	if (simple) {
		int *imap = NULL;
		imap = Calloc(graph->n, int);

		for (i = 0; i < graph->n; i++) {
			imap[i] = i;
		}
		*remap = imap;
		return GMRFLib_SUCCESS;
	}

	/*
	 * check if we have 'global' nodes 
	 */
	for (i = 0, ne = 0; i < graph->n; i++) {
		ne = IMAX(ne, graph->nnbs[i]);
	}
	limit = GMRFLib_GLOBAL_NODE(graph->n, gn_ptr);	       /* this is the limit for a 'global' node */

	if (ne >= limit) {
		/*
		 * yes we have global nodes, make a new graph with these removed. 
		 */
		fixed = Calloc(graph->n, char);
		for (i = 0; i < graph->n; i++) {
			fixed[i] = (graph->nnbs[i] >= limit ? 1 : 0);
		}

		GMRFLib_compute_subgraph(&subgraph, graph, fixed);
		free_subgraph = 1;
	} else {
		subgraph = graph;
		free_subgraph = 0;
	}

	/*
	 * continue with subgraph, which is the original graph minus global nodes 
	 */
	n = subgraph->n;

	if (n > 0) {
		/*
		 * only enter here is the subgraph is non-empty. 
		 */

		for (i = 0, nnz = n; i < n; i++) {
			nnz += subgraph->nnbs[i];
		}

		Q = taucs_ccs_create(n, n, nnz, TAUCS_DOUBLE);
		Q->flags = (TAUCS_PATTERN | TAUCS_SYMMETRIC | TAUCS_TRIANGULAR | TAUCS_LOWER);
		Q->colptr[0] = 0;

		for (i = 0, ic = 0; i < n; i++) {
			Q->rowind[ic++] = i;
			for (k = 0, ne = 1; k < subgraph->nnbs[i]; k++) {
				j = subgraph->nbs[i][k];
				if (j > i) {
					break;
				}
				Q->rowind[ic++] = j;
				ne++;
			}
			Q->colptr[i + 1] = Q->colptr[i] + ne;
		}

		switch (reorder) {
		case GMRFLib_REORDER_IDENTITY:
			p = GMRFLib_strdup("identity");
			break;
		case GMRFLib_REORDER_REVERSE_IDENTITY:
			p = GMRFLib_strdup("reverseidentity");
			break;
		case GMRFLib_REORDER_DEFAULT:
		case GMRFLib_REORDER_METIS:
			p = GMRFLib_strdup("metis");
			break;
		case GMRFLib_REORDER_GENMMD:
			p = GMRFLib_strdup("genmmd");
			break;
		case GMRFLib_REORDER_AMD:
			p = GMRFLib_strdup("amd");
			break;
		case GMRFLib_REORDER_AMDC:
			p = GMRFLib_strdup("amdc");
			break;
		case GMRFLib_REORDER_AMDBAR:
			p = GMRFLib_strdup("amdbar");
			break;
		case GMRFLib_REORDER_AMDBARC:
			p = GMRFLib_strdup("amdbarc");
			break;
		case GMRFLib_REORDER_MD:
			p = GMRFLib_strdup("md");
			break;
		case GMRFLib_REORDER_MMD:
			p = GMRFLib_strdup("mmd");
			break;
		default:
			GMRFLib_ASSERT(0 == 1, GMRFLib_ESNH);
			p = NULL;
		}
		taucs_ccs_order(Q, &perm, &iperm, p);
		Free(p);

		GMRFLib_ASSERT(iperm, GMRFLib_ESNH);
		GMRFLib_ASSERT(perm, GMRFLib_ESNH);

		/*
		 * doit like this to maintain the MEMCHECK facility of GMRFLib 
		 */
		free(perm);
		perm = Calloc(graph->n, int);		       /* yes, need graph->n. */
		memcpy(perm, iperm, n * sizeof(int));
		free(iperm);
		iperm = perm;

		taucs_ccs_free(Q);
	} else {
		/*
		 * in this case, subgraph is empty and we have only global nodes 
		 */
		iperm = NULL;
	}

	if (!free_subgraph) {
		/*
		 * no global nodes, then `iperm' is the reordering 
		 */
		*remap = iperm;				       /* yes, this is correct */
	} else {
		/*
		 * global nodes, correct the reordering computed and add the global nodes at the end 
		 */
		ns = subgraph->n;
		n = graph->n;
		iperm_new = Calloc(graph->n, int);

		for (i = 0; i < ns; i++) {
			iperm_new[subgraph->mothergraph_idx[iperm[i]]] = i;
		}

		/*
		 * in this new code, we sort the global nodes according to the number of neighbours, so the ones with largest number of neighbours are
		 * given highest node-number. 
		 */
		int ng = n - ns;
		int *node = Calloc(ng, int);
		int *nnbs = Calloc(ng, int);

		for (i = 0, j = 0; i < n; i++) {
			if (fixed[i]) {
				node[j] = i;
				nnbs[j] = graph->nnbs[i];
				j++;
			}
		}
		assert(j == ng);

		/*
		 * sort with respect to number of neigbours and carry the node-number along 
		 */
		GMRFLib_qsorts((void *) nnbs, (size_t) ng, sizeof(int), (void *) node, sizeof(int), NULL, 0, GMRFLib_icmp);

		for (i = 0, j = ns; i < ng; i++, j++) {
			iperm_new[node[i]] = j;
		}
		assert(j == n);

		Free(node);
		Free(nnbs);
		GMRFLib_ASSERT(j == n, GMRFLib_ESNH);	       /* just a check... */

		*remap = iperm_new;			       /* this is the reordering */

		GMRFLib_free_graph(subgraph);
		Free(iperm);
	}

	Free(fixed);

	if (!*remap) {
		GMRFLib_ERROR(GMRFLib_EREORDER);
	}

	return GMRFLib_SUCCESS;
}

int GMRFLib_build_sparse_matrix_TAUCS(taucs_ccs_matrix ** L, GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_arg, GMRFLib_graph_tp * graph, int *remap)
{
	int i, j, k, ic, ne, n, nnz, *perm = NULL, *iperm = NULL, id, nan_error = 0;
	taucs_ccs_matrix *Q = NULL;

	id = GMRFLib_thread_id;

	if (!graph || graph->n == 0) {
		*L = NULL;
		return GMRFLib_SUCCESS;
	}

	n = graph->n;
	for (i = 0, nnz = n; i < n; i++) {
		nnz += graph->nnbs[i];
	}

	Q = taucs_ccs_create(n, n, nnz, TAUCS_DOUBLE);
	GMRFLib_ASSERT(Q, GMRFLib_EMEMORY);
	Q->flags = (TAUCS_DOUBLE | TAUCS_SYMMETRIC | TAUCS_TRIANGULAR | TAUCS_LOWER);
	Q->colptr[0] = 0;

#if defined(_OPENMP)
	/*
	 * first, do a first pass to set the indices; then fill the matrix using omp 
	 */
	int *ic_idx = Calloc(n, int);

	for (i = 0, ic = 0; i < n; i++) {
		Q->rowind[ic] = i;
		ic_idx[i] = ic;
		ic++;
		ne = 1;

		for (k = 0; k < graph->nnbs[i]; k++) {
			j = graph->nbs[i][k];
			if (j > i) {
				break;
			}
			Q->rowind[ic] = j;
			ic++;
			ne++;
		}
		Q->colptr[i + 1] = Q->colptr[i] + ne;
	}
#pragma omp parallel for private(i, ic, k, j)
	for (i = 0; i < n; i++) {
		double val;

		ic = ic_idx[i];
		GMRFLib_thread_id = id;

		val = Qfunc(i, i, Qfunc_arg);
		GMRFLib_STOP_IF_NAN_OR_INF(val, i, i);
		Q->values.d[ic++] = val;

		for (k = 0; k < graph->nnbs[i]; k++) {
			j = graph->nbs[i][k];
			if (j > i) {
				break;
			}

			val = Qfunc(i, j, Qfunc_arg);
			GMRFLib_STOP_IF_NAN_OR_INF(val, i, j);
			Q->values.d[ic++] = val;
		}
	}
	GMRFLib_thread_id = id;
	Free(ic_idx);

	if (GMRFLib_catch_error_for_inla) {
		if (nan_error) {
			return !GMRFLib_SUCCESS;
		}
	}
#else
	for (i = 0, ic = 0; i < n; i++) {
		double val;

		Q->rowind[ic] = i;

		val = Qfunc(i, i, Qfunc_arg);
		GMRFLib_STOP_IF_NAN_OR_INF(val, i, i);
		Q->values.d[ic] = val;

		ic++;
		ne = 1;

		for (k = 0; k < graph->nnbs[i]; k++) {
			j = graph->nbs[i][k];
			if (j > i) {
				break;
			}
			Q->rowind[ic] = j;

			val = Qfunc(i, j, Qfunc_arg);
			GMRFLib_STOP_IF_NAN_OR_INF(val, i, j);
			Q->values.d[ic] = val;

			ic++;
			ne++;
		}
		Q->colptr[i + 1] = Q->colptr[i] + ne;
	}
	if (GMRFLib_catch_error_for_inla) {
		if (nan_error) {
			return !GMRFLib_SUCCESS;
		}
	}
#endif
	iperm = remap;					       /* yes, this is correct */
	perm = Calloc(n, int);

	for (i = 0; i < n; i++) {
		perm[iperm[i]] = i;
	}
	*L = taucs_ccs_permute_symmetrically(Q, perm, iperm);  /* permute the matrix */

	taucs_ccs_free(Q);
	Free(perm);

	if (0) {
		static int count = 0;
		char *fnm;
#pragma omp critical
		{
			GMRFLib_sprintf(&fnm, "sparse-matrix-%1d-thread-%1d.txt", count++, omp_get_thread_num());
			FILE *fp = fopen(fnm, "w");
			fprintf(stderr, "write %s\n", fnm);
			for (i = 0; i < n; i++) {
				double qq = Qfunc(i, i, Qfunc_arg);
				fprintf(fp, "%d %d %.16g\n", i, i, qq);
				for (k = 0; k < graph->nnbs[i]; k++) {
					j = graph->nbs[i][k];
					if (i < j) {
						fprintf(fp, "%d %d %.20g\n", i, j, Qfunc(i, j, Qfunc_arg));
					}
				}
			}
			fclose(fp);
			Free(fnm);
		}
	}

	return GMRFLib_SUCCESS;
}

int GMRFLib_factorise_sparse_matrix_TAUCS(taucs_ccs_matrix ** L, supernodal_factor_matrix ** symb_fact, GMRFLib_fact_info_tp * finfo,
					  double **L_inv_diag)
{
	int flags, k, retval;

	if (!L) {
		return GMRFLib_SUCCESS;
	}
	/*
	 * compute some info about the factorization 
	 */
	k = (*L)->colptr[(*L)->n] - (*L)->n;
	finfo->n = (*L)->n;
	finfo->nnzero = 2 * k + (*L)->n;

	flags = (*L)->flags;
	if (!*symb_fact) {
		*symb_fact = (supernodal_factor_matrix *) taucs_ccs_factor_llt_symbolic(*L);
	}

	retval = taucs_ccs_factor_llt_numeric(*L, *symb_fact);
	if (retval) {
		if (GMRFLib_catch_error_for_inla) {
			fprintf(stdout, "\n\t%s\n\tFunction: %s(), Line: %1d, Thread: %1d\n\tFail to factorize Q. I will try to fix it...\n\n",
				RCSId, __GMRFLib_FuncName, __LINE__, omp_get_thread_num());
			return GMRFLib_EPOSDEF;
		} else {
			GMRFLib_ERROR(GMRFLib_EPOSDEF);
		}
	}
	taucs_ccs_free(*L);

	if (include_zeros_in_L) {
		/*
		 * this version will maintain the zero's in L, so that the computation of Qinv gets faster; there is then no need
		 * to check L. 
		 */
		*L = my_taucs_dsupernodal_factor_to_ccs(*symb_fact);
	} else {
		/*
		 * this is the library version which will remove zeros in L. 
		 */
		*L = taucs_supernodal_factor_to_ccs(*symb_fact);
	}
	(*L)->flags = flags & ~TAUCS_SYMMETRIC;		       /* fixes a bug in ver 2.0 av TAUCS */
	taucs_supernodal_factor_free_numeric(*symb_fact);      /* remove the numerics, preserve the symbolic */

	/*
	 * some last info 
	 */
	k = (*L)->colptr[(*L)->n] - (*L)->n;
	finfo->nfillin = k - (finfo->nnzero - finfo->n) / 2;

	/*
	 * compute also the inverse of diag(L) 
	 */
	if (L_inv_diag) {
		int i;

		*L_inv_diag = Calloc((*L)->n, double);
		for (i = 0; i < (*L)->n; i++) {
			(*L_inv_diag)[i] = 1.0 / (*L)->values.d[((*L)->colptr)[i]];
		}
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_free_fact_sparse_matrix_TAUCS(taucs_ccs_matrix * L, double *L_inv_diag, supernodal_factor_matrix * symb_fact)
{
	if (L) {
		taucs_ccs_free(L);
	}
	Free(L_inv_diag);
	if (symb_fact) {
		taucs_supernodal_factor_free(symb_fact);
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_solve_l_sparse_matrix_TAUCS(double *rhs, taucs_ccs_matrix * L, GMRFLib_graph_tp * graph, int *remap)
{
	GMRFLib_EWRAP0(GMRFLib_convert_to_mapped(rhs, NULL, graph, remap));
	GMRFLib_my_taucs_dccs_solve_l(L, rhs);
	GMRFLib_convert_from_mapped(rhs, NULL, graph, remap);
	return GMRFLib_SUCCESS;
}

int GMRFLib_solve_lt_sparse_matrix_TAUCS(double *rhs, taucs_ccs_matrix * L, GMRFLib_graph_tp * graph, int *remap)
{
	double *b = NULL;

	GMRFLib_EWRAP0(GMRFLib_convert_to_mapped(rhs, NULL, graph, remap));
	b = Calloc(graph->n, double);
	memcpy(b, rhs, graph->n * sizeof(double));

	GMRFLib_my_taucs_dccs_solve_lt(L, rhs, b);

	GMRFLib_convert_from_mapped(rhs, NULL, graph, remap);
	Free(b);

	return GMRFLib_SUCCESS;
}

int GMRFLib_solve_llt_sparse_matrix_TAUCS(double *rhs, taucs_ccs_matrix * L, GMRFLib_graph_tp * graph, int *remap)
{
	GMRFLib_EWRAP0(GMRFLib_convert_to_mapped(rhs, NULL, graph, remap));

	if (0) {
		/*
		 * use TAUCS 
		 */
		double *b = Calloc(graph->n, double);
		memcpy(b, rhs, graph->n * sizeof(double));
		taucs_ccs_solve_llt(L, rhs, b);
		Free(b);

		GMRFLib_convert_from_mapped(rhs, NULL, graph, remap);
	} else {
		/*
		 * my version for this particular purpose, a bit faster (15% or so) 
		 */
		GMRFLib_my_taucs_dccs_solve_llt(L, rhs);

		/*
		 * inlined code 
		 */
		int i;
		double *work = Malloc(graph->n, double);

		memcpy(work, rhs, graph->n * sizeof(double));
		for (i = 0; i < graph->n; i++) {
			rhs[i] = work[remap[i]];
		}
		Free(work);
	}

	return GMRFLib_SUCCESS;
}

int GMRFLib_solve_lt_sparse_matrix_special_TAUCS(double *rhs, taucs_ccs_matrix * L, GMRFLib_graph_tp * graph, int *remap, int findx, int toindx,
						 int remapped)
{
	/*
	 * rhs in real world, L in mapped world.  solve L^Tx=b backward only from rhs[findx] up to rhs[toindx].  note that
	 * findx and toindx is in mapped world.  if remapped, do not remap/remap-back the rhs before solving.
	 * 
	 */

	double *b = Calloc(graph->n, double);

	if (!remapped) {
		GMRFLib_convert_to_mapped(rhs, NULL, graph, remap);
	}
	memcpy(&b[toindx], &rhs[toindx], (graph->n - toindx) * sizeof(double));	/* this can be improved */

	GMRFLib_my_taucs_dccs_solve_lt_special(L, rhs, b, findx, toindx);	/* solve it */
	if (!remapped) {
		GMRFLib_convert_from_mapped(rhs, NULL, graph, remap);
	}
	Free(b);

	return GMRFLib_SUCCESS;
}

int GMRFLib_solve_l_sparse_matrix_special_TAUCS(double *rhs, taucs_ccs_matrix * L, GMRFLib_graph_tp * graph, int *remap, int findx, int toindx,
						int remapped)
{
	/*
	 * rhs in real world, L in mapped world.  solve Lx=b backward only from rhs[findx] up to rhs[toindx].  note that
	 * findx and toindx is in mapped world.  if remapped, do not remap/remap-back the rhs before solving.
	 * 
	 */

	double *b = Calloc(graph->n, double);

	if (!remapped) {
		GMRFLib_convert_to_mapped(rhs, NULL, graph, remap);
	}
	memcpy(&b[findx], &rhs[findx], (toindx - findx + 1) * sizeof(double));	/* this can be improved */
	GMRFLib_my_taucs_dccs_solve_l_special(L, rhs, b, findx, toindx);	/* solve it */
	if (!remapped) {
		GMRFLib_convert_from_mapped(rhs, NULL, graph, remap);
	}
	Free(b);

	return GMRFLib_SUCCESS;
}

int GMRFLib_solve_llt_sparse_matrix_special_TAUCS(double *x, taucs_ccs_matrix * L, double *L_inv_diag, GMRFLib_graph_tp * graph, int *remap,
						  int idx)
{
	/*
	 * this is special version of the GMRFLib_solve_llt_sparse_matrix_TAUCS()-routine, where we KNOW that x is 0 exect for a 1 at index
	 * `idx'. return the solution in x. this is requried as this the main task for _ai for large problems. it is easy to switch this
	 * feature off, just modify the GMRFLib_solve_llt_sparse_matrix_special() routine in sparse-interface.c 
	 */

	int n, i, j, ip, jp, idxnew, use_new_code = 1;
	double Aij, Ajj, Aii, *y = NULL, sum;

	GMRFLib_ASSERT(x[idx] == 1.0, GMRFLib_ESNH);

	idxnew = remap[idx];
	x[idx] = 0.0;
	x[idxnew] = 1.0;

	n = L->n;
	y = Calloc(n, double);

	if (use_new_code) {
		GMRFLib_ASSERT(L_inv_diag, GMRFLib_ESNH);
	}

	/*
	 * need only to start at 'idxnew' not 0!; this is the main speedup!!! 
	 */
	if (use_new_code) {
		/*
		 * this version use the L_inv_diag which is 1/diag(L), to simplify some of comptuations, and makes the last expression more compact. 
		 */

		for (j = idxnew; j < n; j++) {
			y[j] = x[j] * L_inv_diag[j];
			for (ip = L->colptr[j] + 1; ip < L->colptr[j + 1]; ip++) {
				i = L->rowind[ip];
				Aij = L->values.d[ip];
				x[i] -= y[j] * Aij;
			}
		}
		for (i = n - 1; i >= 0; i--) {
			sum = 0.0;
			for (jp = L->colptr[i] + 1; jp < L->colptr[i + 1]; jp++) {
				j = L->rowind[jp];
				Aij = L->values.d[jp];
				sum += x[j] * Aij;
			}
			x[i] = (y[i] - sum) * L_inv_diag[i];
		}
	} else {
		for (j = idxnew; j < n; j++) {
			ip = L->colptr[j];
			Ajj = L->values.d[ip];
			y[j] = x[j] / Ajj;

			for (ip = L->colptr[j] + 1; ip < L->colptr[j + 1]; ip++) {
				i = L->rowind[ip];
				Aij = L->values.d[ip];
				x[i] -= y[j] * Aij;
			}
		}
		for (i = n - 1; i >= 0; i--) {
			sum = 0.0;
			for (jp = L->colptr[i] + 1; jp < L->colptr[i + 1]; jp++) {
				j = L->rowind[jp];
				Aij = L->values.d[jp];
				sum += x[j] * Aij;
			}
			y[i] -= sum;
			jp = L->colptr[i];
			Aii = L->values.d[jp];
			x[i] = y[i] / Aii;
		}
	}

	/*
	 * but also reusing y as here, helps 
	 */
	memcpy(y, x, n * sizeof(double));
	for (i = 0; i < n; i++) {
		x[i] = y[remap[i]];
	}
	Free(y);

	return GMRFLib_SUCCESS;
}

int GMRFLib_comp_cond_meansd_TAUCS(double *cmean, double *csd, int indx, double *x, int remapped, taucs_ccs_matrix * L, GMRFLib_graph_tp * graph,
				   int *remap)
{
	/*
	 * compute the conditonal mean and stdev for x[indx]|x[indx+1]...x[n-1] for the current value of x. if `remapped', then 
	 * x is assumed to be remapped for possible (huge) speedup when used repeately.  note: indx is still the user world!!! 
	 */
	int ii = remap[indx];

	if (remapped) {
		GMRFLib_my_taucs_cmsd(cmean, csd, ii, L, x);
	} else {
		GMRFLib_convert_to_mapped(x, NULL, graph, remap);
		GMRFLib_my_taucs_cmsd(cmean, csd, ii, L, x);
		GMRFLib_convert_from_mapped(x, NULL, graph, remap);
	}

	return GMRFLib_SUCCESS;
}

int GMRFLib_log_determinant_TAUCS(double *logdet, taucs_ccs_matrix * L)
{
	int i;

	*logdet = 0.0;
	for (i = 0; i < L->n; i++) {
		*logdet += log(L->values.d[L->colptr[i]]);
	}
	*logdet *= 2;

	return GMRFLib_SUCCESS;
}

int GMRFLib_compute_Qinv_TAUCS(GMRFLib_problem_tp * problem, int storage)
{
	if (!problem) {
		return GMRFLib_SUCCESS;
	}

	if (include_zeros_in_L) {
		/*
		 * we have now modified the code so that zero's maintains in L, so we do not need the checking-step any longer. 
		 */
		GMRFLib_EWRAP0(GMRFLib_compute_Qinv_TAUCS_compute(problem, storage, NULL));
	} else {
		/*
		 * strategy:
		 * 
		 * try first to compute Qinv, if fail, then go back and add zero-terms to L until ok, and compute Qinv again. 
		 */

		int i, n;
		taucs_ccs_matrix *L = NULL, *LL = NULL;	       /* to hold L matrices if new ones are built */
		map_ii **mis_elm = NULL;

		n = problem->sub_sm_fact.TAUCS_L->n;

		/*
		 * no-check, check-once or failsafe? 
		 */
		if ((storage & GMRFLib_QINV_NO_CHECK)) {
			GMRFLib_EWRAP0(GMRFLib_compute_Qinv_TAUCS_compute(problem, storage, NULL));
		} else {
			/*
			 * do some checking 
			 */
			L = GMRFLib_L_duplicate_TAUCS(problem->sub_sm_fact.TAUCS_L, TAUCS_DOUBLE | TAUCS_LOWER | TAUCS_TRIANGULAR);
			while ((mis_elm = GMRFLib_compute_Qinv_TAUCS_check(L))) {
				LL = L;
				L = GMRFLib_compute_Qinv_TAUCS_add_elements(LL, mis_elm);
				taucs_ccs_free(LL);

				for (i = 0; i < n; i++) {
					map_ii_free(mis_elm[i]);
					Free(mis_elm[i]);
				}
				Free(mis_elm);

				if ((storage & GMRFLib_QINV_CHECK_ONCE)) {
					break;		       /* check once only */
				}
			}
			GMRFLib_EWRAP0(GMRFLib_compute_Qinv_TAUCS_compute(problem, storage, L));
			taucs_ccs_free(L);
		}
	}
	return GMRFLib_SUCCESS;
}

map_ii **GMRFLib_compute_Qinv_TAUCS_check(taucs_ccs_matrix * L)
{
	int i, j, k, jp, ii, kk, *nnbs = NULL, **nbs = NULL, *nnbsQ = NULL, nmissing = 0, n;
	map_ii **Qinv_L = NULL;
	map_ii **missing = NULL;

	if (!L) {
		return NULL;
	}
	n = L->n;

	/*
	 * construct a row-list of L_ij's including the diagonal 
	 */
	nnbs = Calloc(n, int);
	nbs = Calloc(n, int *);
	nnbsQ = Calloc(n, int);				       /* number of elm in the Qinv_L[j] hash-table */

	for (i = 0; i < n; i++) {
		for (jp = L->colptr[i]; jp < L->colptr[i + 1]; jp++) {
			nnbs[L->rowind[jp]]++;
		}
	}
#pragma omp parallel for private(i)
	for (i = 0; i < n; i++) {
		nbs[i] = Calloc(nnbs[i], int);
		nnbs[i] = 0;
	}

	for (j = 0; j < n; j++) {
		for (jp = L->colptr[j]; jp < L->colptr[j + 1]; jp++) {	/* including the diagonal */
			i = L->rowind[jp];
			nbs[i][nnbs[i]] = j;
			nnbs[i]++;
			nnbsQ[IMIN(i, j)]++;		       /* for the Qinv_L[] hash-table */
		}
	}

	/*
	 * I join these tasks into one i-loop; preferable for omp 
	 */
	Qinv_L = Calloc(n, map_ii *);
	missing = Calloc(n, map_ii *);
#pragma omp parallel for private(i)
	for (i = 0; i < n; i++) {
		qsort(nbs[i], (size_t) nnbs[i], sizeof(int), GMRFLib_icmp);	/* is this needed ???? */
		Qinv_L[i] = Calloc(1, map_ii);
		map_ii_init_hint(Qinv_L[i], nnbsQ[i]);
		missing[i] = Calloc(1, map_ii);
		map_ii_init(missing[i]);
	}

	if (1) {
		/*
		 * different versions 
		 */
		if (1) {
			/*
			 * fast version 1
			 * 
			 * not much to gain here. 
			 */
			char *Zj = NULL;
			int *Zj_set = NULL, nset;

			Zj = Calloc(n, char);
			Zj_set = Calloc(n, int);

			for (j = n - 1; j >= 0; j--) {
				nset = 0;
				for (k = -1; (k = (int) map_ii_next(Qinv_L[j], k)) != -1;) {
					kk = Qinv_L[j]->contents[k].key;
					Zj[kk] = 1;
					Zj_set[nset++] = kk;
				}

				for (ii = nnbs[j] - 1; ii >= 0; ii--) {
					i = nbs[j][ii];
					for (kk = L->colptr[i] + 1; kk < L->colptr[i + 1]; kk++) {
						k = L->rowind[kk];
						if (Zj[k] == 0 && L->values.d[kk] != 0.0) {
							map_ii_set(missing[IMIN(k, j)], IMAX(k, j), 1);
							if (k < j) {
								map_ii_set(Qinv_L[k], j, 1);	/* also mark those who are missing */
							}
							Zj[k] = 1;
							Zj_set[nset++] = k;
							nmissing++;
						}
					}
					Zj[i] = 1;
					Zj_set[nset++] = i;
					map_ii_set(Qinv_L[i], j, 1);
				}
				/*
				 * if (j > 0) for(kk=0;kk<nset;kk++) Zj[Zj_set[kk]] = 0; 
				 */

				if (j > 0) {		       /* not needed for j=0 */
					if (nset > GMRFLib_NSET_LIMIT(nset, (int) sizeof(char), n)) {
						memset(Zj, 0, n * sizeof(char));
					} else {
						for (kk = 0; kk < nset; kk++) {
							Zj[Zj_set[kk]] = 0;	/* set those to zero */
						}
					}
				}
			}
			Free(Zj);
			Free(Zj_set);
		} else {
			/*
			 * fast version 2
			 * 
			 * run almost as fast as the above, probably since the size of Zj is char and hence small. but is to be 
			 * preferred as it does not matter how large 'nset' is?
			 * 
			 */
			char *Zj = NULL;

			Zj = Calloc(n, char);

			for (j = n - 1; j >= 0; j--) {
				for (k = -1; (k = (int) map_ii_next(Qinv_L[j], k)) != -1;) {
					Zj[Qinv_L[j]->contents[k].key] = 1;
				}

				for (ii = nnbs[j] - 1; ii >= 0; ii--) {
					i = nbs[j][ii];
					for (kk = L->colptr[i] + 1; kk < L->colptr[i + 1]; kk++) {
						k = L->rowind[kk];
						if (Zj[k] == 0 && L->values.d[kk] != 0.0) {
							map_ii_set(missing[IMIN(k, j)], IMAX(k, j), 1);
							if (k < j) {
								map_ii_set(Qinv_L[k], j, 1);	/* also mark those who are missing */
							}
							Zj[k] = 1;
							nmissing++;
						}
					}
					Zj[i] = 1;
					map_ii_set(Qinv_L[i], j, 1);
				}
				if (j > 0) {
					memset(Zj, 0, n * sizeof(char));
				}
			}
			Free(Zj);
		}
	} else {
		/*
		 * slow version, but don't delete! 
		 */
		for (j = n - 1; j >= 0; j--) {
			for (ii = nnbs[j] - 1; ii >= 0; ii--) {
				i = nbs[j][ii];
				for (kk = L->colptr[i] + 1; kk < L->colptr[i + 1]; kk++) {
					k = L->rowind[kk];
					if (!map_ii_ptr(Qinv_L[IMIN(k, j)], IMAX(k, j))) {
						map_ii_set(missing[IMIN(k, j)], IMAX(k, j), 0);
						map_ii_set(Qinv_L[IMIN(k, j)], IMAX(k, j), 0);	/* also mark those who are missing */
						nmissing++;
					}
				}
				map_ii_set(Qinv_L[IMIN(i, j)], IMAX(i, j), 0);
			}
		}
	}

#pragma omp parallel for private(i)
	for (i = 0; i < n; i++) {
		Free(nbs[i]);
		map_ii_free(Qinv_L[i]);
		Free(Qinv_L[i]);
	}
	Free(nbs);
	Free(nnbs);
	Free(Qinv_L);
	Free(nnbsQ);

	if (nmissing == 0) {
		/*
		 * free the missing-hash, and set it to NULL so that this function returns NULL 
		 */
#pragma omp parallel for private(i)
		for (i = 0; i < n; i++) {
			map_ii_free(missing[i]);
			Free(missing[i]);
		}
		Free(missing);
	}

	printf("\n\n\n%s: nmissing %d\n\n", __GMRFLib_FuncName, nmissing);
	if (nmissing) {
		abort();
	}

	return missing;
}

taucs_ccs_matrix *GMRFLib_compute_Qinv_TAUCS_add_elements(taucs_ccs_matrix * L, map_ii ** missing)
{
	int i, j, k, jp, jpp, n, nnz, nnz_new, nmissing, collen;
	taucs_ccs_matrix *LL = NULL;

	n = L->n;
	/*
	 * first count the number of missing terms 
	 */
	for (i = 0, nmissing = 0; i < n; i++) {
		for (k = -1; (k = (int) map_ii_next(missing[i], k)) != -1;) {
			nmissing++;
		}
	}

	nnz = L->colptr[n];
	nnz_new = nnz + nmissing;

	LL = taucs_dccs_create(n, n, nnz_new);
	LL->flags = L->flags;				       /* maintain properties */

	LL->colptr[0] = 0;				       /* start at zero */
	for (j = 0, jp = 0; j < n; j++) {
		/*
		 * do each column j 
		 */
		collen = L->colptr[j + 1] - L->colptr[j];
		jpp = L->colptr[j];

		memcpy(&(LL->rowind[jp]), &(L->rowind[jpp]), collen * sizeof(int));
		memcpy(&(LL->values.d[jp]), &(L->values.d[jpp]), collen * sizeof(double));

		jp += collen;

		/*
		 * now add new elms to this column 
		 */
		for (k = -1; (k = (int) map_ii_next(missing[j], k)) != -1;) {
			LL->rowind[jp] = missing[j]->contents[k].key;
			LL->values.d[jp] = 0.0;
			jp++;
		}
		LL->colptr[j + 1] = jp;
	}

	if (0) {
		for (j = 0; j < n; j++) {
			for (jp = L->colptr[j]; jp < L->colptr[j + 1]; jp++) {
				printf("L[ %1d %1d ] = %.12f\n", L->rowind[jp], j, L->values.d[jp]);
			}
		}

		for (j = 0; j < n; j++) {
			for (jp = LL->colptr[j]; jp < LL->colptr[j + 1]; jp++) {
				printf("LL[ %1d %1d ] = %.12f\n", LL->rowind[jp], j, LL->values.d[jp]);
			}
		}
	}

	return LL;
}
int GMRFLib_compute_Qinv_TAUCS_compute(GMRFLib_problem_tp * problem, int storage, taucs_ccs_matrix * Lmatrix)
{
	/*
	 * compute the elements in Qinv from the non-zero pattern of L (no checking). store them according to `storage':
	 * GMRFLib_QINV_ALL GMRFLib_QINV_NEIGB GMRFLib_QINV_DIAG 
	 */
	double *ptr = NULL, value, diag, *Zj = NULL;
	int i, j, k, jp, ii, kk, jj, iii, jjj, n, *nnbs = NULL, **nbs = NULL, *nnbsQ = NULL, *rremove = NULL, nrremove, *inv_remap =
	    NULL, *Zj_set, nset;
	taucs_ccs_matrix *L = NULL;
	map_ii *mapping = NULL;
	map_id **Qinv_L = NULL, *q = NULL;

	L = (Lmatrix ? Lmatrix : problem->sub_sm_fact.TAUCS_L);	/* chose matrix to use */
	n = L->n;

	/*
	 * construct a row-list of L_ij's including the diagonal 
	 */
	nnbs = Calloc(n, int);
	nbs = Calloc(n, int *);
	nnbsQ = Calloc(n, int);				       /* number of elm in the Qinv_L[j] hash-table */

	for (i = 0; i < n; i++) {
		for (jp = L->colptr[i]; jp < L->colptr[i + 1]; jp++) {
			nnbs[L->rowind[jp]]++;
		}
	}

	for (i = 0; i < n; i++) {
		nbs[i] = Calloc(nnbs[i], int);
		nnbs[i] = 0;
	}

	for (j = 0; j < n; j++) {
		for (jp = L->colptr[j]; jp < L->colptr[j + 1]; jp++) {	/* including the diagonal */
			i = L->rowind[jp];
			nbs[i][nnbs[i]] = j;
			nnbs[i]++;
			nnbsQ[IMIN(i, j)]++;		       /* for the Qinv_L[] hash-table */
		}
	}

	/*
	 * sort and setup the hash-table for storing Qinv_L 
	 */
	Qinv_L = Calloc(n, map_id *);
#pragma omp parallel for private(i)
	for (i = 0; i < n; i++) {
		qsort(nbs[i], (size_t) nnbs[i], sizeof(int), GMRFLib_icmp);	/* needed? */
		Qinv_L[i] = Calloc(1, map_id);
		map_id_init_hint(Qinv_L[i], nnbsQ[i]);
	}

	Zj = Calloc(n, double);
	Zj_set = Calloc(n, int);

	for (j = n - 1; j >= 0; j--) {
		/*
		 * store those indices that are used and set only those to zero 
		 */
		nset = 0;
		q = Qinv_L[j];				       /* just to store the ptr */

		for (k = -1; (k = (int) map_id_next(q, k)) != -1;) {
			jj = q->contents[k].key;
			Zj_set[nset++] = jj;
			Zj[jj] = q->contents[k].value;
		}
		for (ii = nnbs[j] - 1; ii >= 0; ii--) {
			i = nbs[j][ii];
			diag = L->values.d[L->colptr[i]];
			value = (i == j ? 1. / diag : 0.0);
			/*
			 * no gain to omp this loop or to workshare this ii-loop either.... 
			 */
			for (kk = L->colptr[i] + 1; kk < L->colptr[i + 1]; kk++) {
				value -= L->values.d[kk] * Zj[L->rowind[kk]];
			}

			value /= diag;
			Zj[i] = value;
			Zj_set[nset++] = i;

			map_id_set(Qinv_L[i], j, value);
		}
		if (j > 0) {				       /* not needed for j=0 */
			if (nset > GMRFLib_NSET_LIMIT(nset, (int) sizeof(double), n)) {
				memset(Zj, 0, n * sizeof(double));	/* faster if nset is large */
			} else {
				for (kk = 0; kk < nset; kk++) {
					Zj[Zj_set[kk]] = 0.0;  /* set those to zero */
				}
			}
		}
	}
	if (0) {
		/*
		 * keep this OLD version in the source 
		 */
		for (j = n - 1; j >= 0; j--) {
			for (ii = nnbs[j] - 1; ii >= 0; ii--) {
				i = nbs[j][ii];
				diag = L->values.d[L->colptr[i]];
				value = (i == j ? 1. / diag : 0.0);

				for (kk = L->colptr[i] + 1; kk < L->colptr[i + 1]; kk++) {
					k = L->rowind[kk];
					if ((ptr = map_id_ptr(Qinv_L[IMIN(k, j)], IMAX(k, j)))) {
						value -= L->values.d[kk] * *ptr;
					}
				}

				value /= diag;
				map_id_set(Qinv_L[IMIN(i, j)], IMAX(i, j), value);
			}
		}
	}

	/*
	 * compute the mapping 
	 */
	inv_remap = Calloc(n, int);

	for (k = 0; k < n; k++) {
		inv_remap[problem->sub_sm_fact.remap[k]] = k;
	}

	/*
	 * possible remove entries: options are GMRFLib_QINV_ALL GMRFLib_QINV_NEIGB GMRFLib_QINV_DIAG 
	 */
	if (storage & (GMRFLib_QINV_DIAG | GMRFLib_QINV_NEIGB)) {
		rremove = Calloc(n, int);

		for (i = 0; i < n; i++) {
			iii = inv_remap[i];
			if (storage & GMRFLib_QINV_DIAG) {
				for (k = -1, nrremove = 0; (k = (int) map_id_next(Qinv_L[i], k)) != -1;) {
					if ((j = Qinv_L[i]->contents[k].key) != i) {
						rremove[nrremove++] = j;
					}
				}
			} else {
				for (k = -1, nrremove = 0; (k = (int) map_id_next(Qinv_L[i], k)) != -1;) {
					j = Qinv_L[i]->contents[k].key;

					if (j != i) {
						jjj = inv_remap[j];
						if (!GMRFLib_is_neighb(iii, jjj, problem->sub_graph)) {
							rremove[nrremove++] = j;
						}
					}
				}
			}
			for (k = 0; k < nrremove; k++) {
				map_id_remove(Qinv_L[i], rremove[k]);
			}
			map_id_adjustcapacity(Qinv_L[i]);
		}
	}

	/*
	 * correct for constraints, if any. need `iremap' as the matrix terms, constr_m and qi_at_m, is in the sub_graph
	 * coordinates without reordering!
	 * 
	 * not that this is correct for both hard and soft constraints, as the constr_m matrix contains the needed noise-term. 
	 */
	if (problem->sub_constr && problem->sub_constr->nc > 0) {
#pragma omp parallel for private(i, iii, k, j, jjj, kk, value)
		for (i = 0; i < n; i++) {
			iii = inv_remap[i];
			for (k = -1; (k = (int) map_id_next(Qinv_L[i], k)) != -1;) {
				j = Qinv_L[i]->contents[k].key;
				jjj = inv_remap[j];

				map_id_get(Qinv_L[i], j, &value);
				for (kk = 0; kk < problem->sub_constr->nc; kk++) {
					value -= problem->constr_m[iii + kk * n] * problem->qi_at_m[jjj + kk * n];
				}
				map_id_set(Qinv_L[i], j, value);
			}
		}
	}

	/*
	 * done. store Qinv 
	 */
	problem->sub_inverse = Calloc(1, GMRFLib_Qinv_tp);
	problem->sub_inverse->Qinv = Qinv_L;

	/*
	 * compute the mapping for lookup using GMRFLib_Qinv_get(). here, the user lookup using a global index, which is then
	 * transformed to the reordered sub_graph. 
	 */
	problem->sub_inverse->mapping = mapping = Calloc(1, map_ii);
	map_ii_init_hint(mapping, n);
	for (i = 0; i < n; i++) {
		map_ii_set(mapping, problem->sub_graph->mothergraph_idx[i], problem->sub_sm_fact.remap[i]);
	}

	/*
	 * cleanup 
	 */
#pragma omp parallel for private(i)
	for (i = 0; i < n; i++) {
		Free(nbs[i]);
	}
	Free(nbs);
	Free(nnbs);
	Free(nnbsQ);
	Free(inv_remap);
	Free(rremove);
	Free(Zj);
	Free(Zj_set);

	return GMRFLib_SUCCESS;
}

int GMRFLib_my_taucs_dccs_solve_lt(void *vL, double *x, double *b)
{
	taucs_ccs_matrix *L = (taucs_ccs_matrix *) vL;

	int i, j, jp;
	double Aij, Aii;

	for (i = L->n - 1; i >= 0; i--) {
		for (jp = L->colptr[i] + 1; jp < L->colptr[i + 1]; jp++) {
			j = L->rowind[jp];
			Aij = L->values.d[jp];
			b[i] -= x[j] * Aij;
		}

		jp = L->colptr[i];
		Aii = L->values.d[jp];
		x[i] = b[i] / Aii;
	}

	return 0;
}

int GMRFLib_my_taucs_dccs_solve_lt_special(void *vL, double *x, double *b, int from_idx, int to_idx)
{
	taucs_ccs_matrix *L = (taucs_ccs_matrix *) vL;

	int i, j, jp;
	double Aij, Aii;

	for (i = from_idx; i >= to_idx; i--) {
		for (jp = L->colptr[i] + 1; jp < L->colptr[i + 1]; jp++) {
			j = L->rowind[jp];
			Aij = L->values.d[jp];
			b[i] -= x[j] * Aij;
		}

		jp = L->colptr[i];
		Aii = L->values.d[jp];
		x[i] = b[i] / Aii;
	}

	return 0;
}

int GMRFLib_my_taucs_dccs_solve_l_special(void *vL, double *x, double *b, int from_idx, int to_idx)
{
	taucs_ccs_matrix *L = (taucs_ccs_matrix *) vL;
	int ip, i, j;
	double Aij, Ajj;

	for (j = from_idx; j <= to_idx; j++) {
		ip = L->colptr[j];
		Ajj = L->values.d[ip];
		x[j] = b[j] / Ajj;

		for (ip = L->colptr[j] + 1; ip < L->colptr[j + 1]; ip++) {
			i = L->rowind[ip];
			Aij = L->values.d[ip];
			b[i] -= x[j] * Aij;
		}
	}
	return 0;
}

int GMRFLib_my_taucs_dccs_solve_llt(void *vL, double *x)
{
	taucs_ccs_matrix *L = (taucs_ccs_matrix *) vL;
	int n, ip, jp;
	double Aij, Ajj, Aii;

	double *y = NULL;

	n = L->n;

	if (n > 0) {
		y = Calloc(n, double);

		{
			int j;
			for (j = 0; j < n; j++) {
				ip = L->colptr[j];
				Ajj = L->values.d[ip];
				y[j] = x[j] / Ajj;

				for (ip = L->colptr[j] + 1; ip < L->colptr[j + 1]; ip++) {
					int i;

					i = L->rowind[ip];
					Aij = L->values.d[ip];
					x[i] -= y[j] * Aij;
				}
			}
		}

		{
			int i;
			for (i = n - 1; i >= 0; i--) {
				double sum = 0.0;

				for (jp = L->colptr[i] + 1; jp < L->colptr[i + 1]; jp++) {
					int j;

					j = L->rowind[jp];
					Aij = L->values.d[jp];
					sum += x[j] * Aij;
				}
				y[i] -= sum;

				jp = L->colptr[i];
				Aii = L->values.d[jp];
				x[i] = y[i] / Aii;
			}
		}
		Free(y);
	}
	return 0;
}

int GMRFLib_my_taucs_dccs_solve_l(void *vL, double *x)
{
	taucs_ccs_matrix *L = (taucs_ccs_matrix *) vL;
	int n, ip;
	double Aij, Ajj;

	double *y = NULL;

	n = L->n;

	if (n > 0) {
		y = Calloc(n, double);

		{
			int j;
			for (j = 0; j < n; j++) {
				ip = L->colptr[j];
				Ajj = L->values.d[ip];
				y[j] = x[j] / Ajj;

				for (ip = L->colptr[j] + 1; ip < L->colptr[j + 1]; ip++) {
					int i;

					i = L->rowind[ip];
					Aij = L->values.d[ip];
					x[i] -= y[j] * Aij;
				}
			}
		}

		memcpy(x, y, n * sizeof(double));
		Free(y);
	}
	return 0;
}

int GMRFLib_my_taucs_cmsd(double *cmean, double *csd, int idx, taucs_ccs_matrix * L, double *x)
{
	int j, jp;
	double Aij, Aii, b;

	for (b = 0.0, jp = L->colptr[idx] + 1; jp < L->colptr[idx + 1]; jp++) {
		j = L->rowind[jp];
		Aij = L->values.d[jp];
		b -= x[j] * Aij;
	}

	jp = L->colptr[idx];
	Aii = L->values.d[jp];
	*cmean = b / Aii;
	*csd = 1 / Aii;

	return 0;
}

int GMRFLib_my_taucs_check_flags(int flags)
{
#define CheckFLAGS(X) if (flags & X) printf(#X " is ON\n");if (!(flags & X)) printf(#X " is OFF\n")
	CheckFLAGS(TAUCS_INT);
	CheckFLAGS(TAUCS_DOUBLE);
	CheckFLAGS(TAUCS_SINGLE);
	CheckFLAGS(TAUCS_DCOMPLEX);
	CheckFLAGS(TAUCS_SCOMPLEX);
	CheckFLAGS(TAUCS_LOWER);
	CheckFLAGS(TAUCS_UPPER);
	CheckFLAGS(TAUCS_TRIANGULAR);
	CheckFLAGS(TAUCS_SYMMETRIC);
	CheckFLAGS(TAUCS_HERMITIAN);
	CheckFLAGS(TAUCS_PATTERN);
#undef CheckFLAGS
	return GMRFLib_SUCCESS;
}

int GMRFLib_bitmap_factorisation_TAUCS__intern(taucs_ccs_matrix * L, const char *filename)
{
#define ROUND(_i) ((int) ((_i) * reduce_factor))
#define SET(_i, _j) bitmap[ROUND(_i) + ROUND(_j) * N] = 1

	int i, j, jp, n = L->n, N, err;
	double reduce_factor;
	GMRFLib_uchar *bitmap;

	if (GMRFLib_bitmap_max_dimension > 0 && n > GMRFLib_bitmap_max_dimension) {
		N = GMRFLib_bitmap_max_dimension;
		reduce_factor = (double) N / (double) n;
	} else {
		N = n;
		reduce_factor = 1.0;
	}

	bitmap = Calloc(ISQR(N), GMRFLib_uchar);
	for (i = 0; i < n; i++) {
		SET(i, i);
		for (jp = L->colptr[i] + 1; jp < L->colptr[i + 1]; jp++) {
			j = L->rowind[jp];
			SET(i, j);
		}
	}

	err = GMRFLib_bitmap_image(filename, bitmap, N, N);
	Free(bitmap);

	return err;
}

int GMRFLib_bitmap_factorisation_TAUCS(const char *filename_body, taucs_ccs_matrix * L)
{
	/*
	 * create a bitmap-file of the factorization 
	 */
	char *filename = NULL;

	GMRFLib_EWRAP0(GMRFLib_sprintf(&filename, "%s_L.pbm", (filename_body ? filename_body : "taucs_L")));
	GMRFLib_EWRAP0(GMRFLib_bitmap_factorisation_TAUCS__intern(L, filename));
	Free(filename);

	return GMRFLib_SUCCESS;
}

int GMRFLib_amdc(int n, int *pe, int *iw, int *len, int iwlen, int pfree,
		 int *nv, int *next, int *last, int *head, int *elen, int *degree, int ncmpa, int *w)
{
	int result, i;
	double control[AMD_CONTROL], info[AMD_INFO];

	amd_defaults(control);
	result = amd_order(n, pe, iw, last, control, info);
	GMRFLib_ASSERT(result == AMD_OK, GMRFLib_EREORDER);

	for (i = 0; i < n; i++) {
		last[i]++;				       /* to Fortran indexing. */
	}

	return (result == AMD_OK ? GMRFLib_SUCCESS : !GMRFLib_SUCCESS);
}

int GMRFLib_amdbarc(int n, int *pe, int *iw, int *len, int iwlen, int pfree,
		    int *nv, int *next, int *last, int *head, int *elen, int *degree, int ncmpa, int *w)
{
	int result, i;
	double control[AMD_CONTROL], info[AMD_INFO];

	amd_defaults(control);
	control[AMD_AGGRESSIVE] = 0;			       /* turn this off */
	result = amd_order(n, pe, iw, last, control, info);
	GMRFLib_ASSERT(result == AMD_OK, GMRFLib_EREORDER);

	for (i = 0; i < n; i++) {
		last[i]++;				       /* to Fortran indexing. */
	}

	return (result == AMD_OK ? GMRFLib_SUCCESS : !GMRFLib_SUCCESS);
}

#undef GMRFLib_NSET_LIMIT


/* smtp-taucs.c
 * 
 * Copyright (C) 2001-2024 Havard Rue
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

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"
#include "amd.h"
#include "metis.h"

#define TIME_FUNCTIONS 0

GMRFLib_taucs_cache_tp *GMRFLib_taucs_cache_duplicate(GMRFLib_taucs_cache_tp *cache)
{
	if (cache) {
		GMRFLib_taucs_cache_tp *nc = Calloc(1, GMRFLib_taucs_cache_tp);
		nc->n = cache->n;
		nc->nnz = cache->nnz;
		if (nc->n && cache->len) {
			nc->len = Calloc(nc->n, int);
			Memcpy(nc->len, cache->len, nc->n * sizeof(int));
		}
		if (nc->nnz && cache->rowind) {
			nc->rowind = Calloc(nc->nnz, int);
			Memcpy(nc->rowind, cache->rowind, nc->nnz * sizeof(int));
		}
		return nc;
	}
	return NULL;
}

void GMRFLib_taucs_cache_free(GMRFLib_taucs_cache_tp *cache)
{
	if (cache) {
		Free(cache->len);
		Free(cache->rowind);
		Free(cache);
	}
}

taucs_ccs_matrix *my_taucs_dsupernodal_factor_to_ccs(void *vL, GMRFLib_taucs_cache_tp **cache)
{
	GMRFLib_ENTER_ROUTINE;
	supernodal_factor_matrix *L = (supernodal_factor_matrix *) vL;

	int n = L->n;
	if (n == 0) {
		return NULL;
	}

	taucs_ccs_matrix *C = NULL;
	int nnz = 0, *len = NULL;
	int nt = (n < 1E3 ? 1 : (n < 10E4 ? 2 : 4));

	if (cache == NULL || *cache == NULL) {
		len = Calloc(n, int);
		for (int sn = 0; sn < L->n_sn; sn++) {
			int Lsize = L->sn_size[sn];
			int Lup_size = L->sn_up_size[sn];
			int *Lss = L->sn_struct[sn];

			for (int jp = 0; jp < Lsize; jp++) {
				int j = Lss[jp];
				int *len_j = len + j;
				*len_j = 0;

				for (int ip = jp; ip < Lsize; ip++) {
					int i = Lss[ip];
					if (i >= j) {
						(*len_j)++;
						nnz++;
					}
				}

				for (int ip = Lsize; ip < Lup_size; ip++) {
					int i = Lss[ip];
					if (i >= j) {
						(*len_j)++;
						nnz++;
					}
				}
			}
		}
		if (cache) {
			*cache = Calloc(1, GMRFLib_taucs_cache_tp);
			(*cache)->n = n;
			(*cache)->nnz = nnz;
			(*cache)->len = len;
		}
	} else {
		nnz = (*cache)->nnz;
		len = (*cache)->len;
	}

	C = taucs_dccs_create(n, n, nnz);
	C->flags = TAUCS_DOUBLE | TAUCS_TRIANGULAR | TAUCS_LOWER;

	(C->colptr)[0] = 0;
	for (int j = 1; j <= n; j++) {
		(C->colptr)[j] = (C->colptr)[j - 1] + len[j - 1];
	}

	if (cache && (*cache)->rowind) {
		Memcpy(C->rowind, (*cache)->rowind, nnz * sizeof(int));
	} else {
#define CODE_BLOCK							\
		for (int sn = 0; sn < L->n_sn; sn++) {			\
			int *Lss = L->sn_struct[sn];			\
			int Lsize = L->sn_size[sn];			\
			int Lup_size = L->sn_up_size[sn];		\
			for (int jp = 0; jp < Lsize; jp++) {		\
				int j = Lss[jp];			\
				int next = C->colptr[j];		\
				for (int ip = jp; ip < Lsize; ip++) {	\
					int i = Lss[ip];		\
					if (i >= j) {			\
						C->rowind[next++] = i;	\
					}				\
				}					\
				for (int ip = Lsize; ip < Lup_size; ip++) { \
					int i = Lss[ip];		\
					if (i >= j) {			\
						C->rowind[next++] = i;	\
					}				\
				}					\
			}						\
		}

		RUN_CODE_BLOCK(nt, 0, 0);
#undef CODE_BLOCK
		if (cache) {
			(*cache)->rowind = Calloc(nnz, int);
			Memcpy((*cache)->rowind, C->rowind, nnz * sizeof(int));
		}
	}

#define CODE_BLOCK							\
	for (int sn = 0; sn < L->n_sn; sn++) {				\
		int *Lss = L->sn_struct[sn];				\
		int Lsbl = L->sn_blocks_ld[sn];				\
		int Lsize = L->sn_size[sn];				\
		int Lubl = L->up_blocks_ld[sn];				\
		int Lup_size = L->sn_up_size[sn];			\
		taucs_datatype *Lsb = L->sn_blocks[sn];			\
		taucs_datatype *Lub = L->up_blocks[sn];			\
		for (int jp = 0; jp < Lsize; jp++) {			\
			int j = Lss[jp];				\
			int next = C->colptr[j];			\
			taucs_datatype *Lsb_p = Lsb + jp * Lsbl;	\
			taucs_datatype *Lub_p = Lub + jp * Lubl - Lsize; \
			for (int ip = jp; ip < Lsize; ip++) {		\
				int i = Lss[ip];			\
				if (i >= j) {				\
					C->values.d[next++] = Lsb_p[ip]; \
				}					\
			}						\
			for (int ip = Lsize; ip < Lup_size; ip++) {	\
				int i = Lss[ip];			\
				if (i >= j) {				\
					C->values.d[next++] = Lub_p[ip]; \
				}					\
			}						\
		}							\
	}

	RUN_CODE_BLOCK(nt, 0, 0);
#undef CODE_BLOCK

	if (!cache) {
		Free(len);
	}

	GMRFLib_LEAVE_ROUTINE;
	return C;
}

taucs_ccs_matrix *my_taucs_dsupernodal_factor_to_ccs_ORIG(void *vL, GMRFLib_taucs_cache_tp **UNUSED(cache))
{
	GMRFLib_ENTER_ROUTINE;
	// original version, with added unused argument

	supernodal_factor_matrix *L = (supernodal_factor_matrix *) vL;
	taucs_ccs_matrix *C = NULL;
	int n, nnz = 0, *len = NULL;
	taucs_datatype v;

	n = L->n;
	if (n == 0) {
		return NULL;
	}
	len = Calloc(n, int);

	for (int sn = 0; sn < L->n_sn; sn++) {
		int Lsize = L->sn_size[sn];
		int Lup_size = L->sn_up_size[sn];
		int *Lss = L->sn_struct[sn];

		for (int jp = 0; jp < Lsize; jp++) {
			int j = Lss[jp];
			int *len_j = len + j;
			*len_j = 0;
			for (int ip = jp; ip < Lsize; ip++) {
				int i = Lss[ip];
				if (i >= j) {
					(*len_j)++;
					nnz++;
				}
			}
			for (int ip = Lsize; ip < Lup_size; ip++) {
				int i = Lss[ip];
				if (i >= j) {
					(*len_j)++;
					nnz++;
				}
			}
		}
	}

	C = taucs_dccs_create(n, n, nnz);
	C->flags = TAUCS_DOUBLE;
	C->flags |= TAUCS_TRIANGULAR | TAUCS_LOWER;

	(C->colptr)[0] = 0;
	for (int j = 1; j <= n; j++) {
		(C->colptr)[j] = (C->colptr)[j - 1] + len[j - 1];
	}

	for (int sn = 0; sn < L->n_sn; sn++) {
		int *Lss = L->sn_struct[sn];
		int Lsbl = L->sn_blocks_ld[sn];
		int Lsize = L->sn_size[sn];
		int Lubl = L->up_blocks_ld[sn];
		int Lup_size = L->sn_up_size[sn];
		taucs_datatype *Lsb = L->sn_blocks[sn];
		taucs_datatype *Lub = L->up_blocks[sn];

		for (int jp = 0; jp < Lsize; jp++) {
			int j = Lss[jp];
			int next = C->colptr[j];

			taucs_datatype *Lsb_p = Lsb + jp * Lsbl;
			taucs_datatype *Lub_p = Lub + jp * Lubl - Lsize;
			for (int ip = jp; ip < Lsize; ip++) {
				int i = Lss[ip];
				if (i >= j) {
					v = Lsb_p[ip];
					C->rowind[next] = i;
					C->values.d[next] = v;
					next++;
				}
			}
			for (int ip = Lsize; ip < Lup_size; ip++) {
				int i = Lss[ip];
				if (i >= j) {
					v = Lub_p[ip];
					C->rowind[next] = i;
					C->values.d[next] = v;
					next++;
				}
			}
		}
	}

	Free(len);
	GMRFLib_LEAVE_ROUTINE;
	return C;
}

supernodal_factor_matrix *GMRFLib_sm_fact_duplicate_TAUCS(supernodal_factor_matrix *L)
{
#define DUPLICATE(name,len,type) if (1) {				\
		if (L->name && ((len) > 0)) {				\
			LL->name = (type *)Calloc((len), type);		\
			Memcpy(LL->name,L->name,(size_t)(len)*sizeof(type)); \
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

	LL->sn_struct = (int **) Calloc(n_sn, int *);
	for (int i = 0; i < LL->n_sn; i++) {
		DUPLICATE(sn_struct[i], LL->sn_up_size[i], int);
	}

	LL->sn_blocks = (double **) Calloc(n_sn, double *);
	for (int i = 0; i < LL->n_sn; i++) {
		DUPLICATE(sn_blocks[i], ISQR(LL->sn_size[i]), double);
	}

	LL->up_blocks = (double **) Calloc(n_sn, double *);
	for (int i = 0; i < LL->n_sn; i++) {
		DUPLICATE(up_blocks[i], (LL->sn_up_size[i] - LL->sn_size[i]) * (LL->sn_size)[i], double);
	}

#undef DUPLICATE
	return LL;
}

void taucs_ccs_metis5(taucs_ccs_matrix *m, int **perm, int **invperm, char *UNUSED(which))
{
	// this for metis version 5

	int n, nnz;
	int *xadj = NULL;
	int *adj = NULL;
	int *len = NULL;
	int *ptr = NULL;
	int ret;

	if (!(m->flags & TAUCS_SYMMETRIC) && !(m->flags & TAUCS_HERMITIAN)) {
		taucs_printf(Strdup("taucs_ccs_treeorder: METIS ordering only works on symmetric matrices.\n"));
		*perm = NULL;
		*invperm = NULL;
		return;
	}
	/*
	 * this routine may actually work on UPPER as well 
	 */
	if (!(m->flags & TAUCS_LOWER)) {
		taucs_printf(Strdup("taucs_ccs_metis: the lower part of the matrix must be represented.\n"));
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
	for (int i = 0; i < n; i++)
		len[i] = 0;

	for (int j = 0; j < n; j++) {
		for (int ip = (m->colptr)[j]; ip < (m->colptr)[j + 1]; ip++) {
			/*
			 * i = (m->rowind)[ip] - (m->indshift);
			 */
			int i = (m->rowind)[ip];
			if (i != j) {
				len[i]++;
				len[j]++;
			}
		}
	}

	xadj[0] = 0;
	for (int i = 1; i <= n; i++)
		xadj[i] = xadj[i - 1] + len[i - 1];

	for (int i = 0; i < n; i++)
		ptr[i] = xadj[i];

	for (int j = 0; j < n; j++) {
		for (int ip = (m->colptr)[j]; ip < (m->colptr)[j + 1]; ip++) {
			/*
			 * i = (m->rowind)[ip] - (m->indshift);
			 */
			int i = (m->rowind)[ip];
			if (i != j) {
				adj[ptr[i]] = j;
				adj[ptr[j]] = i;
				ptr[i]++;
				ptr[j]++;
			}
		}
	}
	int options[METIS_NOPTIONS];
	// Have to adapt to the PARDISO metis libs
	// METIS_SetDefaultOptions(options);
	for (int i = 0; i < METIS_NOPTIONS; i++) {
		options[i] = -1;
	}

	options[METIS_OPTION_PTYPE] = METIS_PTYPE_KWAY;
	options[METIS_OPTION_CTYPE] = METIS_CTYPE_SHEM;
	options[METIS_OPTION_IPTYPE] = METIS_IPTYPE_NODE;
	options[METIS_OPTION_RTYPE] = METIS_RTYPE_SEP2SIDED;
	options[METIS_OPTION_NCUTS] = 4;
	options[METIS_OPTION_NSEPS] = 4;
	options[METIS_OPTION_NUMBERING] = 0;
	options[METIS_OPTION_SEED] = 2704;
	options[METIS_OPTION_MINCONN] = 1;
	options[METIS_OPTION_NO2HOP] = 0;
	options[METIS_OPTION_COMPRESS] = 1;
	options[METIS_OPTION_PFACTOR] = 200;

	ret = METIS_NodeND(&n, xadj, adj, NULL, options, *perm, *invperm);
	if (ret != METIS_OK) {
		assert(0 == 1);
		return;
	}

	Free(xadj);
	Free(adj);
}

#if defined(INLA_WITH_PARDISO_WORKAROUND)
int METIS51PARDISO_NodeND(int *i, int *j, int *k, int *l, int *m, int *n, int *o)
{
	return METIS_NodeND(i, j, k, l, m, n, o);
}
#endif

size_t GMRFLib_sm_fact_nnz_TAUCS(supernodal_factor_matrix *L)
{
	/*
	 * return the number of non-zeros in the matrix 
	 */
	size_t nnz = 0;
	for (int sn = 0; sn < L->n_sn; sn++) {
		for (int jp = 0; jp < L->sn_size[sn]; jp++) {
			nnz += L->sn_size[sn] - jp;
			nnz += L->sn_up_size[sn] - L->sn_size[sn];
		}
	}
	return (nnz);
}

taucs_ccs_matrix *GMRFLib_L_duplicate_TAUCS(taucs_ccs_matrix *L)
{
	taucs_ccs_matrix *LL = NULL;
	int n, nnz;

	if (!L) {
		return NULL;
	}

	n = L->n;
	nnz = L->colptr[L->n];
	LL = taucs_ccs_create(n, n, nnz, L->flags);

	Memcpy(LL->colptr, L->colptr, (n + 1) * sizeof(int));
	Memcpy(LL->rowind, L->rowind, nnz * sizeof(int));
	Memcpy(LL->values.d, L->values.d, nnz * sizeof(double));

	return LL;
}

taucs_crs_matrix *GMRFLib_LL_duplicate_TAUCS(taucs_crs_matrix *LL)
{
	taucs_crs_matrix *L = NULL;
	int n, nnz;

	if (!LL) {
		return NULL;
	}

	n = LL->n;
	nnz = LL->rowptr[L->n];

	L = Calloc(1, taucs_crs_matrix);
	L->flags = LL->flags;
	L->rowptr = Calloc(n + 1, int);
	L->colind = Calloc(nnz, int);
	L->values.d = Calloc(nnz, double);

	Memcpy(LL->rowptr, L->rowptr, (n + 1) * sizeof(int));
	Memcpy(LL->colind, L->colind, nnz * sizeof(int));
	Memcpy(LL->values.d, L->values.d, nnz * sizeof(double));

	return L;
}

int GMRFLib_print_ccs_matrix(FILE *fp, taucs_ccs_matrix *L)
{
	if (!L) {
		return GMRFLib_SUCCESS;
	}

	int n = L->n;
	int nnz = L->colptr[L->n];

	fprintf(fp, "n = %d\n", n);
	fprintf(fp, "nnz = %d\n", nnz);

	for (int i = 0; i < n + 1; i++) {
		fprintf(fp, "\tcolptr[%1d] = %1d\n", i, L->colptr[i]);
	}
	for (int i = 0; i < nnz; i++) {
		fprintf(fp, "\trowind[%1d] = %1d\n", i, L->rowind[i]);
		fprintf(fp, "\tvalues[%1d] = %.12g\n", i, L->values.d[i]);
	}

	return GMRFLib_SUCCESS;
}

int GMRFLib_compute_reordering_TAUCS(int **remap, GMRFLib_graph_tp *graph, GMRFLib_reorder_tp reorder, GMRFLib_global_node_tp *gn_ptr)
{
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
		if (graph->n >= 0)
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

		GMRFLib_graph_comp_subgraph(&subgraph, graph, fixed, NULL);
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
		Q->flags |= (TAUCS_PATTERN | TAUCS_SYMMETRIC | TAUCS_TRIANGULAR | TAUCS_LOWER);
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
		{
			p = Strdup("identity");
		}
			break;

		case GMRFLib_REORDER_REVERSE_IDENTITY:
		{
			p = Strdup("reverseidentity");
		}
			break;

		case GMRFLib_REORDER_DEFAULT:
		case GMRFLib_REORDER_METIS:
		{
			p = Strdup("metis");
		}
			break;

		case GMRFLib_REORDER_GENMMD:
		{
			p = Strdup("genmmd");
		}
			break;

		case GMRFLib_REORDER_AMD:
		{
			p = Strdup("amd");
		}
			break;

		case GMRFLib_REORDER_AMDC:
		{
			p = Strdup("amdc");
		}
			break;

		case GMRFLib_REORDER_AMDBAR:
		{
			p = Strdup("amdbar");
		}
			break;

		case GMRFLib_REORDER_AMDBARC:
		{
			p = Strdup("amdbarc");
		}
			break;

		case GMRFLib_REORDER_MD:
		{
			p = Strdup("md");
		}
			break;

		case GMRFLib_REORDER_MMD:
		{
			p = Strdup("mmd");
		}
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
		Memcpy(perm, iperm, n * sizeof(int));
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
			iperm_new[iperm[i]] = i;
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
		// GMRFLib_qsort2((void *) nnbs, (size_t) ng, sizeof(int), (void *) node, sizeof(int), NULL, 0, GMRFLib_icmp);
		my_sort2_ii(nnbs, node, ng);

		for (i = 0, j = ns; i < ng; i++, j++) {
			iperm_new[node[i]] = j;
		}
		assert(j == n);

		Free(node);
		Free(nnbs);
		GMRFLib_ASSERT(j == n, GMRFLib_ESNH);	       /* just a check... */

		*remap = iperm_new;			       /* this is the reordering */

		GMRFLib_graph_free(subgraph);
		Free(iperm);
	}

	Free(fixed);

	if (!*remap) {
		GMRFLib_ERROR(GMRFLib_EREORDER);
	}

	return GMRFLib_SUCCESS;
}

int GMRFLib_build_sparse_matrix_TAUCS(int thread_id, taucs_ccs_matrix **L, GMRFLib_Qfunc_tp *Qfunc, void *Qfunc_arg, GMRFLib_graph_tp *graph,
				      int *remap)
{
	int n = 0, *iperm = NULL, nan_error = 0;
	taucs_ccs_matrix *Q = NULL;

	if (!graph || graph->n == 0) {
		*L = NULL;
		return GMRFLib_SUCCESS;
	}

	n = graph->n;
	Q = taucs_ccs_create(n, n, n + graph->nnz, TAUCS_DOUBLE);
	GMRFLib_ASSERT(Q, GMRFLib_EMEMORY);
	Q->flags = (TAUCS_DOUBLE | TAUCS_SYMMETRIC | TAUCS_TRIANGULAR | TAUCS_LOWER);
	Q->colptr[0] = 0;

	GMRFLib_tabulate_Qfunc_arg_tp *arg = (GMRFLib_tabulate_Qfunc_arg_tp *) Qfunc_arg;
	int fast_copy = (Qfunc == GMRFLib_tabulate_Qfunction_std && arg->Q);

	if (fast_copy) {
		Memcpy(Q->rowind, graph->rowidx, (n + graph->nnz / 2) * sizeof(int));
		Memcpy(Q->colptr, graph->colptr, (n + 1) * sizeof(int));
		// for (int i = 0; i < n + graph->nnz / 2; i++) Q->values.d[i] = arg->Q->a[graph->row2col[i]];
		GMRFLib_pack(n + graph->nnz / 2, arg->Q->a, graph->row2col, Q->values.d);
	} else {
		int *ic_idx = Calloc(n, int);
		for (int i = 0, ic = 0; i < n; i++) {
			Q->rowind[ic] = i;
			ic_idx[i] = ic;
			ic++;
			Memcpy(&(Q->rowind[ic]), graph->snbs[i], graph->snnbs[i] * sizeof(int));
			ic += graph->snnbs[i];
			Q->colptr[i + 1] = Q->colptr[i] + graph->snnbs[i] + 1;
		}

#define CODE_BLOCK							\
		for (int i = 0; i < n; i++) {				\
			int ic = ic_idx[i];				\
			double val = Qfunc(thread_id, i, i, NULL, Qfunc_arg);	\
			GMRFLib_STOP_IF_NAN_OR_INF(val, i, i);		\
			Q->values.d[ic++] = val;			\
			for (int k = 0; k < graph->snnbs[i]; k++) {	\
				int j = graph->snbs[i][k];		\
				val = Qfunc(thread_id, i, j, NULL, Qfunc_arg);	\
				GMRFLib_STOP_IF_NAN_OR_INF(val, i, j);	\
				Q->values.d[ic++] = val;		\
			}						\
		}

		RUN_CODE_BLOCK((GMRFLib_Qx_strategy ? GMRFLib_MAX_THREADS() : 1), 0, 0);
#undef CODE_BLOCK

		Free(ic_idx);
	}

	if (nan_error) {
		return !GMRFLib_SUCCESS;
	}

	iperm = remap;					       /* yes, this is correct */

	// 'perm' is not used in the taucs_ccs_permute_symmetrically, so we can just pass NULL

	// int *perm = Calloc(n, int);
	// for (int i = 0; i < n; i++) perm[iperm[i]] = i;

	*L = taucs_ccs_permute_symmetrically(Q, NULL, iperm);  /* permute the matrix */
	taucs_ccs_free(Q);

	return GMRFLib_SUCCESS;
}

int GMRFLib_factorise_sparse_matrix_TAUCS(taucs_ccs_matrix **L, supernodal_factor_matrix **symb_fact, GMRFLib_taucs_cache_tp **cache,
					  GMRFLib_fact_info_tp *finfo)
{
	if (!L) {
		return GMRFLib_SUCCESS;
	}
	assert(*L);

	/*
	 * compute some info about the factorization 
	 */
	int flags, k, retval;
	k = (*L)->colptr[(*L)->n] - (*L)->n;
	finfo->n = (*L)->n;
	finfo->nnzero = 2 * k + (*L)->n;

	flags = (*L)->flags;
	if (!*symb_fact) {
		*symb_fact = (supernodal_factor_matrix *) taucs_ccs_factor_llt_symbolic(*L);
	}

	retval = taucs_ccs_factor_llt_numeric(*L, *symb_fact);
	if (retval) {
		fprintf(stdout, "\n\tFunction: %s(), Line: %1d, Thread: %1d\n\tFailed to factorize Q. I will try to fix it...\n\n",
			__GMRFLib_FuncName, __LINE__, omp_get_thread_num());
		return GMRFLib_EPOSDEF;
	}
	taucs_ccs_free(*L);

	*L = my_taucs_dsupernodal_factor_to_ccs(*symb_fact, cache);
	assert(*L);
	(*L)->flags = flags & ~TAUCS_SYMMETRIC;		       /* fixes a bug in ver 2.0 av TAUCS */
	taucs_supernodal_factor_free_numeric(*symb_fact);      /* remove the numerics, preserve the symbolic */

	/*
	 * some last info 
	 */
	k = (*L)->colptr[(*L)->n] - (*L)->n;
	finfo->nfillin = k - (finfo->nnzero - finfo->n) / 2;

	return GMRFLib_SUCCESS;
}

int GMRFLib_free_fact_sparse_matrix_TAUCS(taucs_ccs_matrix *L, taucs_crs_matrix *LL, supernodal_factor_matrix *symb_fact)
{
	if (L) {
		taucs_ccs_free(L);
	}
	if (LL) {
		taucs_crs_free(LL);
	}
	if (symb_fact) {
		taucs_supernodal_factor_free(symb_fact);
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_solve_l_sparse_matrix_TAUCS(double *rhs, taucs_ccs_matrix *L, GMRFLib_graph_tp *graph, int *remap)
{
	GMRFLib_convert_to_mapped(rhs, NULL, graph, remap);
	GMRFLib_my_taucs_dccs_solve_l(L, rhs);
	GMRFLib_convert_from_mapped(rhs, NULL, graph, remap);
	return GMRFLib_SUCCESS;
}

int GMRFLib_solve_lt_sparse_matrix_TAUCS(double *rhs, taucs_ccs_matrix *L, GMRFLib_graph_tp *graph, int *remap)
{
	static double **wwork = NULL;
	static int *wwork_len = NULL;
	if (!wwork) {
#pragma omp critical (Name_4e65f9abac12404e1d9633582ec69bc86e375bd2)
		{
			if (!wwork) {
				wwork_len = Calloc(GMRFLib_CACHE_LEN(), int);
				wwork = Calloc(GMRFLib_CACHE_LEN(), double *);
			}
		}
	}

	int cache_idx = 0;
	GMRFLib_CACHE_SET_ID(cache_idx);

	if (graph->n > wwork_len[cache_idx]) {
		Free(wwork[cache_idx]);
		wwork_len[cache_idx] = graph->n;
		wwork[cache_idx] = Calloc(wwork_len[cache_idx], double);
	}
	double *work = wwork[cache_idx];
	Memset(work, 0, wwork_len[cache_idx] * sizeof(double));

	double *b = work;
	GMRFLib_convert_to_mapped(rhs, NULL, graph, remap);
	Memcpy(b, rhs, graph->n * sizeof(double));
	GMRFLib_my_taucs_dccs_solve_lt(L, rhs, b);
	GMRFLib_convert_from_mapped(rhs, NULL, graph, remap);

	return GMRFLib_SUCCESS;
}

int GMRFLib_solve_llt_sparse_matrix_TAUCS(double *rhs, taucs_ccs_matrix *L, taucs_crs_matrix *LL, GMRFLib_graph_tp *graph, int *remap, double *work)
{
	assert(graph->n == L->n);
	assert(work);

	GMRFLib_convert_to_mapped(work, rhs, graph, remap);
	if (!LL) {
		GMRFLib_my_taucs_dccs_solve_llt(L, work, rhs);
	} else {
		GMRFLib_my_taucs_dccs_solve_llt3(L, LL, work, rhs);
	}
	GMRFLib_convert_from_mapped(rhs, work, graph, remap);
	return GMRFLib_SUCCESS;
}

int GMRFLib_solve_llt_sparse_matrix2_TAUCS(double *rhs, taucs_ccs_matrix *L, GMRFLib_graph_tp *graph, int *remap, int nrhs, double *work)
{
	// same function but for many rhs.

	int n = graph->n;
	assert(work);

	for (int j = 0; j < nrhs; j++) {
		int offset = j * n;
		GMRFLib_convert_to_mapped(work + offset, rhs + offset, graph, remap);
	}
	GMRFLib_my_taucs_dccs_solve_llt2(L, work, nrhs, rhs);
	for (int j = 0; j < nrhs; j++) {
		int offset = j * n;
		GMRFLib_convert_from_mapped(rhs + offset, work + offset, graph, remap);
	}

	return GMRFLib_SUCCESS;
}

int GMRFLib_solve_lt_sparse_matrix_special_TAUCS(double *rhs, taucs_ccs_matrix *L, GMRFLib_graph_tp *graph, int *remap, int findx, int toindx,
						 int remapped)
{
	/*
	 * rhs in real world, L in mapped world.  solve L^Tx=b backward only from rhs[findx] up to rhs[toindx].  note that
	 * findx and toindx is in mapped world.  if remapped, do not remap/remap-back the rhs before solving.
	 * 
	 */

	static double **wwork = NULL;
	static int *wwork_len = NULL;
	if (!wwork) {
#pragma omp critical (Name_9c6d559b5470558ef474f5640951d6b63990a46d)
		{
			if (!wwork) {
				wwork_len = Calloc(GMRFLib_CACHE_LEN(), int);
				wwork = Calloc(GMRFLib_CACHE_LEN(), double *);
			}
		}
	}

	int cache_idx = 0;
	GMRFLib_CACHE_SET_ID(cache_idx);

	if (graph->n > wwork_len[cache_idx]) {
		Free(wwork[cache_idx]);
		wwork_len[cache_idx] = graph->n;
		wwork[cache_idx] = Calloc(wwork_len[cache_idx], double);
	}
	double *work = wwork[cache_idx];
	Memset(work, 0, wwork_len[cache_idx] * sizeof(double));

	double *b = work;
	if (!remapped) {
		GMRFLib_convert_to_mapped(rhs, NULL, graph, remap);
	}
	Memcpy(&b[toindx], &rhs[toindx], (graph->n - toindx) * sizeof(double));	/* this can be improved */

	GMRFLib_my_taucs_dccs_solve_lt_special(L, rhs, b, findx, toindx);	/* solve it */
	if (!remapped) {
		GMRFLib_convert_from_mapped(rhs, NULL, graph, remap);
	}

	return GMRFLib_SUCCESS;
}

int GMRFLib_solve_l_sparse_matrix_special_TAUCS(double *rhs, taucs_ccs_matrix *L, GMRFLib_graph_tp *graph, int *remap, int findx, int toindx,
						int remapped)
{
	/*
	 * rhs in real world, L in mapped world.  solve Lx=b backward only from rhs[findx] up to rhs[toindx].  note that
	 * findx and toindx is in mapped world.  if remapped, do not remap/remap-back the rhs before solving.
	 */

	static double **wwork = NULL;
	static int *wwork_len = NULL;
	if (!wwork) {
#pragma omp critical (Name_a3dba7d9a29b2dbf1981362774e31bd1c94148ec)
		{
			if (!wwork) {
				wwork_len = Calloc(GMRFLib_CACHE_LEN(), int);
				wwork = Calloc(GMRFLib_CACHE_LEN(), double *);
			}
		}
	}

	int cache_idx = 0;
	GMRFLib_CACHE_SET_ID(cache_idx);

	if (graph->n > wwork_len[cache_idx]) {
		Free(wwork[cache_idx]);
		wwork_len[cache_idx] = graph->n;
		wwork[cache_idx] = Calloc(wwork_len[cache_idx], double);
	}
	double *work = wwork[cache_idx];
	Memset(work, 0, wwork_len[cache_idx] * sizeof(double));

	double *b = work;
	if (!remapped) {
		GMRFLib_convert_to_mapped(rhs, NULL, graph, remap);
	}
	Memcpy(&b[findx], &rhs[findx], (toindx - findx + 1) * sizeof(double));	/* this can be improved */
	GMRFLib_my_taucs_dccs_solve_l_special(L, rhs, b, findx, toindx);	/* solve it */
	if (!remapped) {
		GMRFLib_convert_from_mapped(rhs, NULL, graph, remap);
	}

	return GMRFLib_SUCCESS;
}

int GMRFLib_solve_llt_sparse_matrix_special_TAUCS(double *x, taucs_ccs_matrix *L, GMRFLib_graph_tp *UNUSED(graph), int *remap, int idx)
{
	/*
	 * this is special version of the GMRFLib_solve_llt_sparse_matrix_TAUCS()-routine, where we KNOW that x is 0 exect for a 1 at index
	 * `idx'. return the solution in x. this is requried as this the main task for _ai for large problems.
	 */

	static double **wwork = NULL;
	static int *wwork_len = NULL;
	if (!wwork) {
#pragma omp critical (Name_ae25603ba826d85ac7ffa0b88a9f11d5c2246a83)
		{
			if (!wwork) {
				wwork_len = Calloc(GMRFLib_CACHE_LEN(), int);
				wwork = Calloc(GMRFLib_CACHE_LEN(), double *);
			}
		}
	}

	int cache_idx = 0;
	GMRFLib_CACHE_SET_ID(cache_idx);
	GMRFLib_ASSERT(x[idx] == 1.0, GMRFLib_ESNH);

	int idx_new = remap[idx];
	int n = L->n;
	x[idx] = 0.0;
	x[idx_new] = 1.0;

	if (n > wwork_len[cache_idx]) {
		Free(wwork[cache_idx]);
		wwork_len[cache_idx] = n;
		wwork[cache_idx] = Calloc(wwork_len[cache_idx], double);
	}
	double *work = wwork[cache_idx];
	double *y = work;
	GMRFLib_fill(idx_new, 0.0, y);

	double *d = L->values.d;
	int *colptr = L->colptr;
	int *rowind = L->rowind;

	for (int j = idx_new; j < n; j++) {
		y[j] = x[j] / d[colptr[j]];
		double yj = -y[j];
		for (int ip = colptr[j] + 1; ip < colptr[j + 1]; ip++) {
			int i = rowind[ip];
			x[i] = fma(yj, d[ip], x[i]);
		}
	}
	for (int i = n - 1; i >= 0; i--) {
		int jp = colptr[i];
		int jp1 = colptr[i] + 1;
		double inv_Aii = 1.0 / d[jp];
		double s = GMRFLib_ddot_idx_mkl(colptr[i + 1] - jp1, d + jp1, x, rowind + jp1);
		x[i] = (y[i] - s) * inv_Aii;
	}

	Memcpy(y, x, n * sizeof(double));
	GMRFLib_pack(n, y, remap, x);

	return GMRFLib_SUCCESS;
}

int GMRFLib_comp_cond_meansd_TAUCS(double *cmean, double *csd, int indx, double *x, int remapped, taucs_ccs_matrix *L, GMRFLib_graph_tp *graph,
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

int GMRFLib_log_determinant_TAUCS(double *logdet, taucs_ccs_matrix *L)
{
	*logdet = 0.0;
	for (int i = 0; i < L->n; i++) {
		*logdet += log(L->values.d[L->colptr[i]]);
	}
	*logdet *= 2.0;

	return GMRFLib_SUCCESS;
}

int GMRFLib_compute_Qinv_TAUCS(GMRFLib_problem_tp *problem)
{
	if (!problem) {
		return GMRFLib_SUCCESS;
	}

	GMRFLib_compute_Qinv_TAUCS_compute(problem, NULL);
	return GMRFLib_SUCCESS;
}

int GMRFLib_compute_Qinv_TAUCS_compute(GMRFLib_problem_tp *problem, taucs_ccs_matrix *Lmatrix)
{
	int n = 0, *nnbs = NULL, **nbs = NULL, *nnbsQ = NULL, *inv_remap = NULL;
	taucs_ccs_matrix *L = NULL;
	map_id **Qinv_L = NULL;

	L = (Lmatrix ? Lmatrix : problem->sub_sm_fact.TAUCS_L);	/* chose matrix to use */
	n = L->n;
	assert(n >= 0);

	/*
	 * construct a row-list of L_ij's including the diagonal 
	 */

	iCalloc_init(2 * n, 2);
	nnbs = iCalloc_get(n);
	nnbsQ = iCalloc_get(n);

	for (int i = 0; i < n; i++) {
		for (int jp = L->colptr[i]; jp < L->colptr[i + 1]; jp++) {
			nnbs[L->rowind[jp]]++;
		}
	}

	int mm = GMRFLib_isum(n, nnbs);
	int *work_nnbs = Calloc(mm, int);
	nbs = Calloc(n, int *);
	for (int i = 0, m = 0; i < n; m += nnbs[i], i++) {
		nbs[i] = work_nnbs + m;
	}
	Memset(nnbs, 0, n * sizeof(int));

	for (int j = 0; j < n; j++) {
		for (int jp = L->colptr[j]; jp < L->colptr[j + 1]; jp++) {	/* including the diagonal */
			int i = L->rowind[jp];
			nbs[i][nnbs[i]] = j;
			nnbs[i]++;
			nnbsQ[IMIN(i, j)]++;		       /* for the Qinv_L[] hash-table */
		}
	}

	/*
	 * sort and setup the hash-table for storing Qinv_L 
	 */
	Qinv_L = Calloc(n, map_id *);
#pragma omp parallel for
	for (int i = 0; i < n; i++) {
		GMRFLib_qsort(nbs[i], (size_t) nnbs[i], sizeof(int), GMRFLib_icmp);
		Qinv_L[i] = Calloc(1, map_id);
		map_id_init_hint(Qinv_L[i], nnbsQ[i]);
	}

	double *Zj = Calloc(n, double);
	double *d = L->values.d;
	for (int j = n - 1; j >= 0; j--) {
		// store those indices that are used and set only those to zero 
		map_id *q = Qinv_L[j];
		for (int k = -1; (k = (int) map_id_next(q, k)) != -1;) {
			int jj = q->contents[k].key;
			Zj[jj] = q->contents[k].value;
		}

		for (int ii = nnbs[j] - 1; ii >= 0; ii--) {
			int i = nbs[j][ii];
			int nn = L->colptr[i + 1] - (L->colptr[i] + 1);
			int kk = L->colptr[i] + 1;
			double diag = L->values.d[L->colptr[i]];
			double dot = GMRFLib_ddot_idx_mkl(nn, d + kk, Zj, L->rowind + kk);
			double value = (i == j ? 1.0 / diag : 0.0);
			value = (value - dot) / diag;
			Zj[i] = value;
			map_id_set(Qinv_L[i], j, value);
		}
	}

	// compute the mapping 
	inv_remap = Calloc(n, int);
	for (int k = 0; k < n; k++) {
		inv_remap[problem->sub_sm_fact.remap[k]] = k;
	}

	// its good to remove as then we do not need to correct that many for constraints
	int *rremove = nnbsQ;
	Memset(rremove, 0, n * sizeof(int));
	for (int i = 0; i < n; i++) {
		int iii = inv_remap[i];
		int nrremove = 0;
		for (int k = -1; (k = (int) map_id_next(Qinv_L[i], k)) != -1;) {
			int j = Qinv_L[i]->contents[k].key;
			if (j != i) {
				int jjj = inv_remap[j];
				if (!GMRFLib_graph_is_nb(iii, jjj, problem->sub_graph)) {
					rremove[nrremove++] = j;
				}
			}
		}
		for (int k = 0; k < nrremove; k++) {
			map_id_remove(Qinv_L[i], rremove[k]);
		}
		// this can be costly, so we ignore. this will do 'realloc'
		// map_id_adjustcapacity(Qinv_L[i]);
	}

	/*
	 * correct for constraints, if any. need `iremap' as the matrix terms, constr_m and qi_at_m, is in the sub_graph
	 * coordinates without reordering!
	 * 
	 * not that this is correct for both hard and soft constraints, as the constr_m matrix contains the needed noise-term. 
	 */
	if (problem->sub_constr && problem->sub_constr->nc > 0) {
#pragma omp parallel for
		for (int i = 0; i < n; i++) {
			int inc = n;
			int iii = inv_remap[i];
			double *xx = &(problem->constr_m[iii]);
			for (int k = -1; (k = (int) map_id_next(Qinv_L[i], k)) != -1;) {
				int j = Qinv_L[i]->contents[k].key;
				int jjj = inv_remap[j];
				double value = 0.0;
				double *yy = &(problem->qi_at_m[jjj]);
				double sum = ddot_(&(problem->sub_constr->nc), xx, &inc, yy, &inc);
				map_id_get(Qinv_L[i], j, &value);
				map_id_set(Qinv_L[i], j, value - sum);
			}
		}
	}

	problem->sub_inverse = Calloc(1, GMRFLib_Qinv_tp);
	problem->sub_inverse->Qinv = Qinv_L;

	/*
	 * compute the mapping for lookup using GMRFLib_Qinv_get(). here, the user lookup using a global index, which is then
	 * transformed to the reordered sub_graph. 
	 */
	problem->sub_inverse->mapping = Calloc(n, int);
	Memcpy(problem->sub_inverse->mapping, problem->sub_sm_fact.remap, n * sizeof(int));

	Free(inv_remap);
	Free(Zj);
	Free(work_nnbs);
	Free(nbs);
	iCalloc_free();

	return GMRFLib_SUCCESS;
}

int GMRFLib_my_taucs_dccs_solve_lt(void *vL, double *x, double *b)
{
	taucs_ccs_matrix *L = (taucs_ccs_matrix *) vL;

	for (int i = L->n - 1; i >= 0; i--) {
		int jp1 = L->colptr[i] + 1;
		b[i] -= GMRFLib_ddot_idx_mkl(L->colptr[i + 1] - jp1, L->values.d + jp1, x, L->rowind + jp1);

		int jp = L->colptr[i];
		double Aii = L->values.d[jp];
		x[i] = b[i] / Aii;
	}

	return 0;
}

int GMRFLib_my_taucs_dccs_solve_lt_special(void *vL, double *x, double *b, int from_idx, int to_idx)
{
	taucs_ccs_matrix *L = (taucs_ccs_matrix *) vL;

	GMRFLib_fill(from_idx, 0.0, x);
	for (int i = from_idx; i >= to_idx; i--) {
		int jp1 = L->colptr[i] + 1;
		b[i] -= GMRFLib_ddot_idx_mkl(L->colptr[i + 1] - jp1, L->values.d + jp1, x, L->rowind + jp1);

		int jp = L->colptr[i];
		double Aii = L->values.d[jp];
		x[i] = b[i] / Aii;
	}

	return 0;
}

int GMRFLib_my_taucs_dccs_solve_l_special(void *vL, double *x, double *b, int from_idx, int to_idx)
{
	taucs_ccs_matrix *L = (taucs_ccs_matrix *) vL;
	for (int j = from_idx; j <= to_idx; j++) {
		int ip = L->colptr[j];
		double Ajj = L->values.d[ip];
		x[j] = b[j] / Ajj;
		for (ip = L->colptr[j] + 1; ip < L->colptr[j + 1]; ip++) {
			int i = L->rowind[ip];
			double Aij = L->values.d[ip];
			b[i] -= x[j] * Aij;
		}
	}
	return 0;
}

int GMRFLib_my_taucs_dccs_solve_llt(void *__restrict vL, double *__restrict x, double *__restrict w)
{
	taucs_ccs_matrix *L = (taucs_ccs_matrix *) vL;
	int n = L->n;
	if (n == 0) {
		return 0;
	}
#if TIME_FUNCTIONS
	static double tref[] = { 0, 0 };
	static double count = 0;
#endif

	assert(w);
	double *y = w;

	double *d = L->values.d;
	int *colptr = L->colptr;
	int *rowind = L->rowind;

#if TIME_FUNCTIONS
	tref[0] -= GMRFLib_timer();
#endif

	int jfirst = GMRFLib_find_nonzero(x, n, 1);
	if (jfirst < 0) {
		// only zero's
		GMRFLib_fill(n, 0.0, x);
		return 0;
	} else {
		GMRFLib_fill(jfirst, 0.0, y);
	}

	for (int j = jfirst; j < n; j++) {
		y[j] = x[j] / d[colptr[j]];
		double yj = -y[j];
		for (int ip = colptr[j] + 1; ip < colptr[j + 1]; ip++) {
			int i = rowind[ip];
			x[i] = fma(yj, d[ip], x[i]);
		}
	}

#if TIME_FUNCTIONS
	tref[0] += GMRFLib_timer();
	tref[1] -= GMRFLib_timer();
#endif

	for (int i = n - 1; i >= 0; i--) {
		int jp = colptr[i];
		int jp1 = jp + 1;
		double inv_Aii = 1.0 / d[jp];
		y[i] -= GMRFLib_ddot_idx_mkl(colptr[i + 1] - jp1, d + jp1, x, rowind + jp1);
		x[i] = y[i] * inv_Aii;
	}

#if TIME_FUNCTIONS
	tref[1] += GMRFLib_timer();
	count++;
	printf("solve_llt:  relative %.6g %.6g abs %g %g %g\n", tref[0] / (tref[0] + tref[1]), tref[1] / (tref[0] + tref[1]),
	       tref[0] / count, tref[1] / count, (tref[0] + tref[1]) / count);
#endif

	return 0;
}

int GMRFLib_my_taucs_dccs_solve_llt2(void *__restrict vL, double *__restrict x, int nrhs, double *__restrict w)
{
	taucs_ccs_matrix *L = (taucs_ccs_matrix *) vL;
	int n = L->n;

	if (n <= 0 || nrhs <= 0) {
		return 0;
	}
#if TIME_FUNCTIONS
	static double tref[] = { 0, 0 };
	static double count = 0;
#endif

	double *work = w;
	int ione = 1;

#if TIME_FUNCTIONS
	tref[0] -= GMRFLib_timer();
#endif

	// check the case where the rhs contains 0's from the beginning. then we can start at the first non-zero index
	int jfirst = n;					       /* do not change, see below */
	for (int k = 0; k < nrhs; k++) {
		int jj = GMRFLib_find_nonzero(x + k * n, n, 1);
		if (jj >= 0) {
			jfirst = IMIN(jfirst, jj);
		}
		if (ISZERO(jfirst)) {
			break;
		}
	}

	if (jfirst == n) {
		GMRFLib_fill(n * nrhs, 0.0, x);
		return 0;
	}

	Memcpy(work, x, n * nrhs * sizeof(double));
	for (int j = 0; j < nrhs; j++) {
		double *xx = x + j;
		double *ww = work + j * n;
		dcopy_(&n, ww, &ione, xx, &nrhs);
	}

	double *y = work;
	GMRFLib_fill(nrhs * jfirst, 0.0, y);
	for (int j = jfirst; j < n; j++) {
		int ip = L->colptr[j];
		int offset_j = j * nrhs;
		double iAjj = 1.0 / L->values.d[ip];
		double *yy = y + offset_j;
		double *xx = x + offset_j;

#pragma omp simd
		for (int k = 0; k < nrhs; k++) {
			yy[k] = xx[k] * iAjj;
		}

		for (ip = L->colptr[j] + 1; ip < L->colptr[j + 1]; ip++) {
			double Aij = -L->values.d[ip];	       // OOOPS! add minus here for daxpy
			int offset_i = L->rowind[ip] * nrhs;
			xx = x + offset_i;
			GMRFLib_daxpy(nrhs, Aij, yy, xx);
		}
	}

#if TIME_FUNCTIONS
	tref[0] += GMRFLib_timer();
	tref[1] -= GMRFLib_timer();
#endif

	for (int i = n - 1; i >= 0; i--) {
		double sum[nrhs];
		GMRFLib_fill(nrhs, 0.0, sum);
		for (int jp = L->colptr[i] + 1; jp < L->colptr[i + 1]; jp++) {
			int offset_j = L->rowind[jp] * nrhs;
			double Aij = L->values.d[jp];
			double *xx = x + offset_j;
			GMRFLib_daxpy(nrhs, Aij, xx, sum);
		}

		int offset_i = i * nrhs;
		double *yy = y + offset_i;
		GMRFLib_daxpy(nrhs, -1.0, sum, yy);

		int jp = L->colptr[i];
		double iAii = 1.0 / L->values.d[jp];
		double *xx = x + offset_i;
		yy = y + offset_i;
#pragma omp simd
		for (int k = 0; k < nrhs; k++) {
			xx[k] = yy[k] * iAii;
		}
	}

	Memcpy(work, x, n * nrhs * sizeof(double));
	for (int j = 0; j < nrhs; j++) {
		double *xx = x + j * n;
		double *ww = work + j;
		dcopy_(&n, ww, &nrhs, xx, &ione);
	}

#if TIME_FUNCTIONS
	tref[1] += GMRFLib_timer();
	count += nrhs;
	printf("solve_llt2:  relative %.6g %.6g abs %g %g %g\n", tref[0] / (tref[0] + tref[1]), tref[1] / (tref[0] + tref[1]),
	       tref[0] / count, tref[1] / count, (tref[0] + tref[1]) / count);
#endif

	return 0;
}

int GMRFLib_my_taucs_dccs_solve_llt3(void *vL, void *vLL, double *x, double *w)
{
	// this version using both ccs and crs. 

	taucs_ccs_matrix *L = (taucs_ccs_matrix *) vL;
	taucs_crs_matrix *LL = (taucs_crs_matrix *) vLL;
	int n = L->n;

#if TIME_FUNCTIONS
	static double tref[] = { 0, 0 };
	static double count = 0;
#endif

	if (n == 0) {
		return 0;
	}

	assert(w);
	double *y = w;
	double *d = LL->values.d;
	int *rowptr = LL->rowptr;
	int *colind = LL->colind;

#if TIME_FUNCTIONS
	tref[0] -= GMRFLib_timer();
#endif

	y[0] = x[0] / d[0];
	for (int i = 1; i < n; i++) {
		int m = rowptr[i + 1] - rowptr[i];
		int jj = rowptr[i];
		double s = GMRFLib_ddot_idx_mkl(m, d + jj, y, colind + jj);
		y[i] = (x[i] - s) / d[rowptr[i + 1] - 1];
	}

#if TIME_FUNCTIONS
	tref[0] += GMRFLib_timer();
	tref[1] -= GMRFLib_timer();
#endif

	d = L->values.d;
	int *colptr = L->colptr;
	int *rowind = L->rowind;

	for (int i = n - 1; i >= 0; i--) {
		int jp = colptr[i];
		int jp1 = jp + 1;
		double inv_Aii = 1.0 / d[jp];
		y[i] -= GMRFLib_ddot_idx_mkl(colptr[i + 1] - jp1, d + jp1, x, rowind + jp1);
		x[i] = y[i] * inv_Aii;
	}

#if TIME_FUNCTIONS
	tref[1] += GMRFLib_timer();
	count++;
	printf("solve_llt3:  relative %.6g %.6g abs %g %g %g\n", tref[0] / (tref[0] + tref[1]), tref[1] / (tref[0] + tref[1]),
	       tref[0] / count, tref[1] / count, (tref[0] + tref[1]) / count);
#endif

	return 0;
}

int GMRFLib_my_taucs_dccs_solve_l(void *vL, double *x)
{
	taucs_ccs_matrix *L = (taucs_ccs_matrix *) vL;
	int n = L->n;

	static double **wwork = NULL;
	static int *wwork_len = NULL;
	if (!wwork) {
#pragma omp critical (Name_adb454feb2a421a0a2effd2a5298f308a1c3f192)
		{
			if (!wwork) {
				wwork_len = Calloc(GMRFLib_CACHE_LEN(), int);
				wwork = Calloc(GMRFLib_CACHE_LEN(), double *);
			}
		}
	}

	int cache_idx = 0;
	GMRFLib_CACHE_SET_ID(cache_idx);

	if (n > wwork_len[cache_idx]) {
		Free(wwork[cache_idx]);
		wwork_len[cache_idx] = n;
		wwork[cache_idx] = Calloc(wwork_len[cache_idx], double);
	}
	double *work = wwork[cache_idx];
	Memset(work, 0, wwork_len[cache_idx] * sizeof(double));

	double *y = work;
	if (n > 0) {
		for (int j = 0; j < n; j++) {
			int ip = L->colptr[j];
			double Ajj = L->values.d[ip];
			y[j] = x[j] / Ajj;
			double yj = -y[j];
			for (ip = L->colptr[j] + 1; ip < L->colptr[j + 1]; ip++) {
				double Aij = L->values.d[ip];
				int i = L->rowind[ip];
				x[i] = fma(yj, Aij, x[i]);
			}
		}
		Memcpy(x, y, n * sizeof(double));
	}
	return 0;
}

int GMRFLib_my_taucs_cmsd(double *cmean, double *csd, int idx, taucs_ccs_matrix *L, double *x)
{
	double b = 0.0;
	for (int jp = L->colptr[idx] + 1; jp < L->colptr[idx + 1]; jp++) {
		int j = L->rowind[jp];
		double Aij = L->values.d[jp];
		b -= x[j] * Aij;
	}

	int jp = L->colptr[idx];
	double Aii = L->values.d[jp];
	*cmean = b / Aii;
	*csd = 1.0 / Aii;

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

int GMRFLib_bitmap_factorisation_TAUCS__intern(taucs_ccs_matrix *L, const char *filename)
{
#define ROUND(_i) ((int) ((_i) * reduce_factor))
#define SET(_i, _j) bitmap[ROUND(_i) + ROUND(_j) * N] = 1

	int n = L->n, N, err;
	double reduce_factor;
	GMRFLib_uchar *bitmap = NULL;

	if (GMRFLib_bitmap_max_dimension > 0 && n > GMRFLib_bitmap_max_dimension) {
		N = GMRFLib_bitmap_max_dimension;
		reduce_factor = (double) N / (double) n;
	} else {
		N = n;
		reduce_factor = 1.0;
	}

	bitmap = Calloc(ISQR(N), GMRFLib_uchar);
	for (int i = 0; i < n; i++) {
		SET(i, i);
		for (int jp = L->colptr[i] + 1; jp < L->colptr[i + 1]; jp++) {
			int j = L->rowind[jp];
			SET(i, j);
		}
	}

	err = GMRFLib_bitmap_image(filename, bitmap, N, N);
	Free(bitmap);

	return err;
}

int GMRFLib_bitmap_factorisation_TAUCS(const char *filename_body, taucs_ccs_matrix *L)
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

int GMRFLib_amdc(int n, int *pe, int *iw, int *UNUSED(len), int UNUSED(iwlen), int UNUSED(pfree),
		 int *UNUSED(nv), int *UNUSED(next), int *last, int *UNUSED(head), int *UNUSED(elen),
		 int *UNUSED(degree), int UNUSED(ncmpa), int *UNUSED(w))
{
	int result;
	double control[AMD_CONTROL], info[AMD_INFO];

	amd_defaults(control);
	result = amd_order(n, pe, iw, last, control, info);
	GMRFLib_ASSERT(result == AMD_OK, GMRFLib_EREORDER);

	for (int i = 0; i < n; i++) {
		last[i]++;				       /* to Fortran indexing. */
	}

	return (result == AMD_OK ? GMRFLib_SUCCESS : !GMRFLib_SUCCESS);
}

int GMRFLib_amdbarc(int n, int *pe, int *iw, int *UNUSED(len), int UNUSED(iwlen), int UNUSED(pfree),
		    int *UNUSED(nv), int *UNUSED(next), int *last, int *UNUSED(head), int *UNUSED(elen),
		    int *UNUSED(degree), int UNUSED(ncmpa), int *UNUSED(w))
{
	int result;
	double control[AMD_CONTROL], info[AMD_INFO];

	amd_defaults(control);
	control[AMD_AGGRESSIVE] = 0;			       /* turn this off */
	result = amd_order(n, pe, iw, last, control, info);
	GMRFLib_ASSERT(result == AMD_OK, GMRFLib_EREORDER);

	for (int i = 0; i < n; i++) {
		last[i]++;				       /* to Fortran indexing. */
	}

	return (result == AMD_OK ? GMRFLib_SUCCESS : !GMRFLib_SUCCESS);
}

taucs_crs_matrix *GMRFLib_ccs2crs(taucs_ccs_matrix *L)
{
	GMRFLib_ENTER_ROUTINE;

	const int debug = 0;
	taucs_crs_matrix *LL = Calloc(1, taucs_crs_matrix);

	LL->n = L->n;
	LL->m = L->m;
	LL->flags = L->flags;

	int n = L->n;
	int nnz = L->colptr[n];

	LL->rowptr = Calloc(n + 1, int);
	LL->colind = Calloc(nnz, int);
	LL->values.d = Calloc(nnz, double);

	// number of elements pr column
	int *clen = Calloc(n, int);

	for (int j = 0; j < n; j++) {
		int ip = L->colptr[j];
		clen[j]++;
		for (ip = L->colptr[j] + 1; ip < L->colptr[j + 1]; ip++) {
			int i = L->rowind[ip];
			clen[i]++;
		}
	}

	for (int i = 1; i <= n; i++) {
		LL->rowptr[i] = LL->rowptr[i - 1] + clen[i - 1];
	}

	// reuse storage with a different name
	int *rowidx = clen;
	Memset(rowidx, 0, n * sizeof(int));

	for (int j = 0; j < n; j++) {
		int ip = L->colptr[j];
		double Ajj = L->values.d[ip];

		int k = LL->rowptr[j + 1] - 1;
		LL->colind[k] = j;
		LL->values.d[k] = Ajj;

		for (ip = L->colptr[j] + 1; ip < L->colptr[j + 1]; ip++) {
			int i = L->rowind[ip];
			double Aij = L->values.d[ip];
			k = LL->rowptr[i] + rowidx[i];
			LL->colind[k] = j;
			LL->values.d[k] = Aij;
			rowidx[i]++;
		}
	}

#define CODE_BLOCK							\
	for (int i = 0; i < n; i++) {					\
		int m = LL->rowptr[i + 1] - LL->rowptr[i];		\
		int j = LL->rowptr[i];					\
		my_sort2_id(LL->colind + j, (double *) LL->values.d + j, m); \
	}

	RUN_CODE_BLOCK(IMIN(4, GMRFLib_MAX_THREADS()), 0, 0);
#undef CODE_BLOCK

	if (debug) {
		printf("CCS\n");
		printf("colptr ");
		for (int i = 0; i < n + 1; i++) {
			printf(" %1d", L->colptr[i]);
		}
		printf("\nrowind ");
		for (int i = 0; i < nnz; i++) {
			printf(" %1d", L->rowind[i]);
		}
		printf("\nvalues ");
		for (int i = 0; i < nnz; i++) {
			printf(" %.2f", L->values.d[i]);
		}
		printf("\n");

		printf("CRS\n");
		printf("rowptr ");
		for (int i = 0; i < n + 1; i++) {
			printf(" %1d", LL->rowptr[i]);
		}
		printf("\ncolind ");
		for (int i = 0; i < nnz; i++) {
			printf(" %1d", LL->colind[i]);
		}
		printf("\nvalues ");
		for (int i = 0; i < nnz; i++) {
			printf(" %.2f", LL->values.d[i]);
		}
		printf("\n");
	}

	Free(rowidx);
	GMRFLib_LEAVE_ROUTINE;
	return (LL);
}

void taucs_crs_free(taucs_crs_matrix *L)
{
	Free(L->rowptr);
	Free(L->colind);
	Free(L->values.d);
	Free(L);
}

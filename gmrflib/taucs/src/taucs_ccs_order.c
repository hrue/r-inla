#ifdef TAUCS_CORE_GENERAL
#       include <limits.h>
#       if __has_include(<malloc.h>)
#              include <malloc.h>
#       endif
#       include <stdio.h>
#       include <stdlib.h>
#       include <string.h>
#       include <assert.h>
#       include <math.h>
#       include "taucs.h"
#       include "../external/colamd.h"
int GMRFLib_amdc(int n, int *pe, int *iw, int *len, int iwlen, int pfree, int *nv, int *next, int *last, int *head, int *elen, int *degree,
		 int ncmpa, int *w);
int GMRFLib_amdbarc(int n, int *pe, int *iw, int *len, int iwlen, int pfree, int *nv, int *next, int *last, int *head, int *elen, int *degree,
		    int ncmpa, int *w);
static void taucs_ccs_colamd(taucs_ccs_matrix *m, int **perm, int **invperm, char *UNUSED(which))
{
#       ifndef TAUCS_CONFIG_COLAMD
	*perm = NULL;
	*invperm = NULL;
	return;
#       else
	double knobs[COLAMD_KNOBS];
	int Alen = 0;
	int *A = NULL;
	int *p = NULL;
	int *ip = NULL;
	int k, nnz;
	int i;
	if (m->flags & TAUCS_SYMMETRIC || m->flags & TAUCS_HERMITIAN) {
		return;
	}
	*perm = NULL;
	*invperm = NULL;
	nnz = (m->colptr)[m->n];
	p = (int *) taucs_malloc((m->n + 1) * sizeof(int));
	ip = (int *) taucs_malloc((m->n + 1) * sizeof(int));
	assert(p && ip);
	Alen = colamd_recommended(nnz, m->m, m->n);
	if (Alen >= 0)
		A = taucs_malloc(Alen * sizeof(int));
	assert(A);
	assert(A);
	colamd_set_defaults(knobs);
	for (i = 0; i <= m->n; i++)
		p[i] = (m->colptr)[i];
	for (k = 0; k < nnz; k++)
		A[k] = (m->rowind)[k];
	if (!colamd(m->m, m->n, Alen, A, p, knobs)) {
		taucs_free(A);
		taucs_free(p);
		return;
	}
	taucs_free(A);
	*perm = p;
	*invperm = ip;
	for (i = 0; i < m->n; i++)
		(*invperm)[(*perm)[i]] = i;
#       endif
}
static void taucs_ccs_amd(taucs_ccs_matrix *m, int **perm, int **invperm, char *which)
{
#       ifndef TAUCS_CONFIG_AMD
	*perm = NULL;
	*invperm = NULL;
	return;
#       else
	int n = 0, iwlen = 0, pfree = 0, ncmpa = 0, iovflo = 0;
	int *iw = NULL;
	int *pe = NULL;
	int *degree = NULL;
	int *nv = NULL;
	int *next = NULL;
	int *last = NULL;
	int *head = NULL;
	int *elen = NULL;
	int *w = NULL;
	int *len = NULL;
	int nnz = 0, i = 0, j = 0, ip = 0;
	if (!(m->flags & TAUCS_SYMMETRIC) && !(m->flags & TAUCS_HERMITIAN)) {
		*perm = NULL;
		*invperm = NULL;
		return;
	}
	if (!(m->flags & TAUCS_LOWER)) {
		*perm = NULL;
		*invperm = NULL;
		return;
	}
	*perm = NULL;
	*invperm = NULL;
	n = m->n;
	nnz = (m->colptr)[n];
	pe = (int *) taucs_malloc((n + 1) * sizeof(int));
	degree = (int *) taucs_malloc(n * sizeof(int));
	nv = (int *) taucs_malloc(n * sizeof(int));
	next = (int *) taucs_malloc(n * sizeof(int));
	last = (int *) taucs_malloc(n * sizeof(int));
	head = (int *) taucs_malloc(n * sizeof(int));
	elen = (int *) taucs_malloc(n * sizeof(int));
	w = (int *) taucs_malloc(n * sizeof(int));
	len = (int *) taucs_malloc(n * sizeof(int));
	iwlen = n + (int) (10.0 * (nnz - n));
	iw = (int *) taucs_malloc(iwlen * sizeof(int));
	if (!pe || !degree || !nv || !next || !last || !head || !elen || !w || !len || !iw) {
		taucs_free(pe);
		taucs_free(degree);
		taucs_free(nv);
		taucs_free(next);
		taucs_free(last);
		taucs_free(head);
		taucs_free(elen);
		taucs_free(w);
		taucs_free(len);
		taucs_free(iw);
		return;
	}
	int offset;
	if (!strcmp(which, "amdc") || !strcmp(which, "amdbarc")) {
		offset = 0;				       /* C */
	} else {
		offset = 1;				       /* Fortran */
	}
	iovflo = INT_MAX;				       /* for 32-bit only! */
	for (i = 0; i < n; i++)
		len[i] = 0;
	for (j = 0; j < n; j++) {
		for (ip = (m->colptr)[j]; ip < (m->colptr)[j + 1]; ip++) {
			i = (m->rowind)[ip];
			if (i != j) {
				len[i]++;
				len[j]++;
			}
		}
	}
	pe[0] = offset;
	for (i = 1; i < n; i++)
		pe[i] = pe[i - 1] + len[i - 1];
	pfree = pe[n - 1] + len[n - 1];
	if (offset == 0) {
		pe[n] = pfree;
	}
	/*
	 * use degree as a temporary 
	 */
	for (i = 0; i < n; i++)
		degree[i] = pe[i] - offset;
	for (j = 0; j < n; j++) {
		for (ip = (m->colptr)[j]; ip < (m->colptr)[j + 1]; ip++) {
			i = (m->rowind)[ip];
			if (i != j) {
				iw[degree[i]] = j + offset;
				iw[degree[j]] = i + offset;
				degree[i]++;
				degree[j]++;
			}
		}
	}
	if (!strcmp(which, "mmd"))
		amdexa_(&n, pe, iw, len, &iwlen, &pfree, nv, next, last, head, elen, degree, &ncmpa, w, &iovflo);
	else if (!strcmp(which, "md"))
		amdtru_(&n, pe, iw, len, &iwlen, &pfree, nv, next, last, head, elen, degree, &ncmpa, w, &iovflo);
	else if (!strcmp(which, "amdbar"))
		amdbarnew_(&n, pe, iw, len, &iwlen, &pfree, nv, next, last, head, elen, degree, &ncmpa, w);
	else if (!strcmp(which, "amd"))
		amdnew_(&n, pe, iw, len, &iwlen, &pfree, nv, next, last, head, elen, degree, &ncmpa, w);
	else if (!strcmp(which, "amdc"))
		GMRFLib_amdc(n, pe, iw, len, iwlen, pfree, nv, next, last, head, elen, degree, ncmpa, w);
	else if (!strcmp(which, "amdbarc"))
		GMRFLib_amdbarc(n, pe, iw, len, iwlen, pfree, nv, next, last, head, elen, degree, ncmpa, w);
	else {
		return;
	}
	taucs_free(pe);
	taucs_free(degree);
	taucs_free(nv);
	taucs_free(next);
	taucs_free(head);
	taucs_free(elen);
	taucs_free(w);
	taucs_free(iw);
	for (i = 0; i < n; i++)
		last[i]--;
	for (i = 0; i < n; i++)
		len[last[i]] = i;
	*perm = last;
	*invperm = len;
#       endif
}
static void taucs_ccs_genmmd(taucs_ccs_matrix *m, int **perm, int **invperm, char *UNUSED(which))
{
#       ifndef TAUCS_CONFIG_GENMMD
	*perm = NULL;
	*invperm = NULL;
	return;
#       else
	int n, maxint, delta, nofsub;
	int *xadj;
	int *adjncy;
	int *invp;
	int *prm;
	int *dhead;
	int *qsize;
	int *llist;
	int *marker;
	int *len;
	int *next;
	int nnz, i, j, ip;
	if (!(m->flags & TAUCS_SYMMETRIC) && !(m->flags & TAUCS_HERMITIAN)) {
		*perm = NULL;
		*invperm = NULL;
		return;
	}
	if (!(m->flags & TAUCS_LOWER)) {
		*perm = NULL;
		*invperm = NULL;
		return;
	}
	*perm = NULL;
	*invperm = NULL;
	n = m->n;
	nnz = (m->colptr)[n];
	delta = 1;					       /* DELTA is a parameter to allow the choice of nodes whose degree <= min-degree +
							        * DELTA. */
	delta = 1;					       /* DELTA is a parameter to allow the choice of nodes whose degree <= min-degree +
							        * DELTA. */
	maxint = 32000;
	assert(sizeof(int) == 4);
	maxint = 2147483647;				       /* 2**31-1, for 32-bit only! */
	xadj = (int *) taucs_malloc((n + 1) * sizeof(int));
	adjncy = (int *) taucs_malloc((2 * nnz - n) * sizeof(int));
	invp = (int *) taucs_malloc((n + 1) * sizeof(int));
	prm = (int *) taucs_malloc(n * sizeof(int));
	dhead = (int *) taucs_malloc((n + 1) * sizeof(int));
	qsize = (int *) taucs_malloc((n + 1) * sizeof(int));
	llist = (int *) taucs_malloc(n * sizeof(int));
	marker = (int *) taucs_malloc(n * sizeof(int));
	if (!xadj || !adjncy || !invp || !prm || !dhead || !qsize || !llist || !marker) {
		taucs_free(xadj);
		taucs_free(adjncy);
		taucs_free(invp);
		taucs_free(prm);
		taucs_free(dhead);
		taucs_free(qsize);
		taucs_free(llist);
		taucs_free(marker);
		return;
	}
	len = dhead;					       /* we reuse space */
	next = qsize;					       /* we reuse space */
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
			} else {
				/*
				 * len[i] ++;
				 */
			}
		}
	}
	xadj[0] = 1;
	for (i = 1; i <= n; i++)
		xadj[i] = xadj[i - 1] + len[i - 1];
	for (i = 0; i < n; i++)
		next[i] = xadj[i] - 1;
	for (j = 0; j < n; j++) {
		for (ip = (m->colptr)[j]; ip < (m->colptr)[j + 1]; ip++) {
			/*
			 * i = (m->rowind)[ip] - (m->indshift);
			 */
			i = (m->rowind)[ip];
			assert(next[i] < 2 * nnz - n);
			assert(next[j] < 2 * nnz - n);
			if (i != j) {
				adjncy[next[i]] = j + 1;
				adjncy[next[j]] = i + 1;
				next[i]++;
				next[j]++;
			} else {
				/*
				 * adjncy[ next[i] ] = j+1; next[i] ++; 
				 */
			}
		}
	}
	genmmd_(&n, xadj, adjncy, invp, prm, &delta, dhead, qsize, llist, marker, &maxint, &nofsub);
	taucs_free(marker);
	taucs_free(llist);
	taucs_free(qsize);
	taucs_free(dhead);
	taucs_free(xadj);
	taucs_free(adjncy);
	for (i = 0; i < n; i++)
		prm[i]--;
	for (i = 0; i < n; i++)
		invp[prm[i]] = i;
	*perm = prm;
	*invperm = invp;
#       endif
}
static void taucs_ccs_treeorder(taucs_ccs_matrix *m, int **perm, int **invperm)
{
	int n, nnz, i, j, ip, k, p, nleaves;
	int *adjptr;
	int *adj;
	int *len;
	int *ptr;
	int *degree;
	int *leaves;
	if (!(m->flags & TAUCS_SYMMETRIC) && !(m->flags & TAUCS_HERMITIAN)) {
		*perm = NULL;
		*invperm = NULL;
		return;
	}
	/*
	 * this routine may actually work on UPPER as well 
	 */
	if (!(m->flags & TAUCS_LOWER)) {
		*perm = NULL;
		*invperm = NULL;
		return;
	}
	n = m->n;
	nnz = (m->colptr)[n];
	*perm = (int *) taucs_malloc(n * sizeof(int));
	*invperm = (int *) taucs_malloc(n * sizeof(int));
	len = (int *) taucs_malloc(n * sizeof(int));
	degree = (int *) taucs_malloc(n * sizeof(int));
	leaves = (int *) taucs_malloc(n * sizeof(int));
	adjptr = (int *) taucs_malloc(n * sizeof(int));
	adj = (int *) taucs_malloc(2 * (nnz - n) * sizeof(int));
	if (!(*perm) || !(*invperm) || !adjptr || !adj || !len || !degree || !leaves) {
		taucs_free(adj);
		taucs_free(adjptr);
		taucs_free(len);
		taucs_free(leaves);
		taucs_free(degree);
		taucs_free(*perm);
		taucs_free(*invperm);
		*perm = *invperm = NULL;
	}
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
	nleaves = 0;
	for (i = 0; i < n; i++) {
		degree[i] = len[i];
		if (degree[i] <= 1) {
			leaves[nleaves] = i;
			nleaves++;
		}
	}
	adjptr[0] = 0;
	for (i = 1; i < n; i++)
		adjptr[i] = adjptr[i - 1] + len[i - 1];
	ptr = *perm;
	for (i = 0; i < n; i++)
		ptr[i] = adjptr[i];
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
	for (i = 0; i < n; i++) {
		nleaves--;
		if (nleaves <= 0) {
			taucs_free(adj);
			taucs_free(adjptr);
			taucs_free(len);
			taucs_free(leaves);
			taucs_free(degree);
			taucs_free(*perm);
			taucs_free(*invperm);
			*perm = *invperm = NULL;
		}
		j = leaves[nleaves];
		(*perm)[i] = j;
		(*invperm)[j] = i;
		if (len[j] > 0) {
			if (len[j] != 1) {
				/*
				 * not a tree 
				 */
				taucs_free(adj);
				taucs_free(adjptr);
				taucs_free(len);
				taucs_free(leaves);
				taucs_free(degree);
				taucs_free(*perm);
				taucs_free(*invperm);
				*perm = *invperm = NULL;
			}
			p = adj[adjptr[j]];
			for (k = 0; k < len[p]; k++)
				if (adj[adjptr[p] + k] == j)
					break;
			if (k >= len[p]) {		       /* otherwise j does not show up in p's adjacency list */
				taucs_free(adj);
				taucs_free(adjptr);
				taucs_free(len);
				taucs_free(leaves);
				taucs_free(degree);
				taucs_free(*perm);
				taucs_free(*invperm);
				*perm = *invperm = NULL;
			}
			len[p]--;
			for (; k < len[p]; k++)
				adj[adjptr[p] + k] = adj[adjptr[p] + k + 1];
			if (len[p] == 1) {		       /* degree was higher and now is 1 */
				leaves[nleaves] = p;
				nleaves++;
			}
		}
	}
	taucs_free(adj);
	taucs_free(adjptr);
	taucs_free(len);
	taucs_free(leaves);
	taucs_free(degree);
}
void METIS_NodeND(int *, int *, int *, int *, int *, int *, int *);
void METIS51PARDISO_NodeND(int *, int *, int *, int *, int *, int *, int *);
#       if !defined(USE_METIS4)
void taucs_ccs_metis5(taucs_ccs_matrix * m, int **perm, int **invperm, char *which);
#       endif
static void taucs_ccs_metis(taucs_ccs_matrix *m, int **perm, int **invperm, char *UNUSED(which))
{
	// this for metis version 4
#       ifndef TAUCS_CONFIG_METIS
	*perm = NULL;
	*invperm = NULL;
	return;
#       else
	int n, nnz, i, j, ip;
	int *xadj;
	int *adj;
	int num_flag = 0;
	int options_flag[8];
	int *len;
	int *ptr;
	if (!(m->flags & TAUCS_SYMMETRIC) && !(m->flags & TAUCS_HERMITIAN)) {
		*perm = NULL;
		*invperm = NULL;
		return;
	}
	if (!(m->flags & TAUCS_LOWER)) {
		*perm = NULL;
		*invperm = NULL;
		return;
	}
	n = m->n;
	nnz = (m->colptr)[n];
	*perm = (int *) taucs_malloc(n * sizeof(int));
	*invperm = (int *) taucs_malloc(n * sizeof(int));
	xadj = (int *) taucs_malloc((n + 1) * sizeof(int));
	adj = (int *) taucs_malloc(2 * nnz * sizeof(int));
	if (!(*perm) || !(*invperm) || !xadj || !adj) {
		taucs_free(*perm);
		taucs_free(*invperm);
		taucs_free(xadj);
		taucs_free(adj);
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
	options_flag[0] = 1;				       /* use these options */
	options_flag[1] = 3;				       /* default */
	options_flag[2] = 1;				       /* default */
	options_flag[3] = 1;				       /* two-side refinement */
	options_flag[4] = 0;				       /* no debug */
	options_flag[5] = 1;				       /* default */
	options_flag[6] = 0;				       /* this is slow if non-zero. global nodes */
	options_flag[7] = 3;				       /* number of separators */
#              if defined(NO_PARDISO_LIB)
	METIS_NodeND(&n, xadj, adj, &num_flag, options_flag, *perm, *invperm);
#              else
	METIS51PARDISO_NodeND(&n, xadj, adj, &num_flag, options_flag, *perm, *invperm);
#              endif
	taucs_free(xadj);
	taucs_free(adj);
#       endif
}
static void taucs_ccs_randomperm(int n, int **perm, int **invperm)
{
	int i;
	*perm = (int *) taucs_malloc(n * sizeof(int));
	*invperm = (int *) taucs_malloc(n * sizeof(int));
	if (!(*perm) || !(*invperm)) {
		taucs_free(*perm);
		taucs_free(*invperm);
		*perm = *invperm = NULL;
		return;
	}
	for (i = 0; i < n; i++)
		(*perm)[i] = i;
	for (i = 0; i < n; i++) {
		int i1, i2;
		int t;
		i1 = rand() % (n - i);
		i2 = n - i - 1;
		t = (*perm)[i1];
		(*perm)[i1] = (*perm)[i2];
		(*perm)[i2] = t;
	}
	for (i = 0; i < n; i++)
		(*invperm)[(*perm)[i]] = i;
	return;
}
void taucs_ccs_order(taucs_ccs_matrix *m, int **perm, int **invperm, char *which)
{
	if (!strcmp(which, "mmd") || !strcmp(which, "amd") || !strcmp(which, "md") || !strcmp(which, "amdbar") ||
	    !strcmp(which, "amdc") || !strcmp(which, "amdbarc"))
		taucs_ccs_amd(m, perm, invperm, which);
	else if (!strcmp(which, "metis")) {
#       if defined(USE_METIS4)
		taucs_ccs_metis(m, perm, invperm, which);
#       else
		taucs_ccs_metis5(m, perm, invperm, which);
#       endif
	} else if (!strcmp(which, "genmmd"))
		taucs_ccs_genmmd(m, perm, invperm, which);
	else if (!strcmp(which, "colamd"))
		taucs_ccs_colamd(m, perm, invperm, which);
	else if (!strcmp(which, "random"))
		taucs_ccs_randomperm(m->n, perm, invperm);
	else if (!strcmp(which, "tree")) {
		taucs_ccs_treeorder(m, perm, invperm);
		if (*perm == NULL)			       /* perhaps the graph of the matrix is not a tree */
			taucs_ccs_metis(m, perm, invperm, "metis");
	} else if (!strcmp(which, "identity")) {
		int i;
		*perm = (int *) taucs_malloc((m->n) * sizeof(int));
		*invperm = (int *) taucs_malloc((m->n) * sizeof(int));
		if (!(*perm) || !(*invperm)) {
			taucs_free(*perm);
			taucs_free(*invperm);
			*perm = *invperm = NULL;
			return;
		}
		for (i = 0; i < m->n; i++)
			(*perm)[i] = (*invperm)[i] = i;
		return;
	} else {
		*perm = *invperm = NULL;
	}
}
#endif							       /* TAUCS_CORE_GENERAL */

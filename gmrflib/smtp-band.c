
/* GMRFLib-smtp-band.c
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
  \file smtp-band.c
  \brief The implementation of the interface towards LAPACK's band solver.
*/

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

#ifndef HGVERSION
#define HGVERSION
#endif
static const char RCSId[] = "file: " __FILE__ "  " HGVERSION;

/* Pre-hg-Id: $Id: smtp-band.c,v 1.47 2010/02/26 17:55:22 hrue Exp $ */

int GMRFLib_compute_reordering_BAND(int **remap, GMRFLib_graph_tp * graph)
{
	/*
	 * compute the reordering from the graph using the routine in acm582.F 
	 */
	int i, j, lconnec, bandwidth, profile, error, space, ioptpro, worklen, *rstart, *connec, *degree, *work, simple;

	if (!graph || !graph->n)
		return GMRFLib_SUCCESS;

	/*
	 * check if we have a simple solution --> no neigbours 
	 */
	for (i = 0, simple = 1; i < graph->n && simple; i++) {
		simple = (graph->nnbs[i] > 0 ? 0 : 1);
	}
	if (simple) {
		int *imap;
		imap = Calloc(graph->n, int);

		for (i = 0; i < graph->n; i++) {
			imap[i] = i;
		}
		*remap = imap;
		return GMRFLib_SUCCESS;
	}

	/*
	 * task I, reformat the graph to fit the fortran-routines 
	 */
	lconnec = 0;
	for (i = 0; i < graph->n; i++) {
		lconnec += graph->nnbs[i];
	}
	if (lconnec == 0) {
		/*
		 * no connections in the graph, use the identity-map. 
		 */
		*remap = Calloc(graph->n, int);

		for (i = 0; i < graph->n; i++) {
			(*remap)[i] = i;
		}
		return GMRFLib_SUCCESS;
	}

	connec = Calloc(lconnec, int);
	rstart = Calloc(graph->n, int);

	degree = graph->nnbs;				       /* yes! */

	rstart[0] = 1;					       /* fortran indx'ing */
	for (i = 1; i < graph->n; i++) {
		rstart[i] = degree[i - 1] + rstart[i - 1];
	}
	for (i = 0; i < graph->n; i++) {
		for (j = 0; j < degree[i]; j++) {
			connec[rstart[i] - 1 + j] = graph->nbs[i][j] + 1;	/* fortran indx'ing */
		}
	}

	worklen = 6 * graph->n + 3;			       /* maximum over all graphs */
	work = Calloc(worklen, int);

	*remap = Calloc(graph->n, int);

	for (i = 0; i < graph->n; i++) {
		(*remap)[i] = i + 1;			       /* fortran indx'ing */
	}

	ioptpro = 0;					       /* do bandwidth reduction */
	error = space = 0;
	gpskca_(&graph->n, degree, rstart, connec, &ioptpro, &worklen, *remap, work, &bandwidth, &profile, &error, &space);

	for (i = 0; i < graph->n; i++) {
		(*remap)[i]--;				       /* correct for fortran indx'ing */
	}

	Free(work);
	Free(rstart);
	Free(connec);

	if (error) {
		GMRFLib_ERROR(GMRFLib_EREORDER);
	}

	return GMRFLib_SUCCESS;
}

int GMRFLib_build_sparse_matrix_BAND(double **bandmatrix, GMRFLib_Qfunc_tp * Qfunc, void *Qfunc_arg, GMRFLib_graph_tp * graph, int *remap,
				     int bandwidth)
{
#define BIDX(i,j) ((i)+(j)*nrow)			       /* band index'ing */

	/*
	 * return a band-matrix BMATRIX in L-storage defining the precision matrix 
	 */

	int i, ncol, nrow, id, nan_error = 0;

	id = GMRFLib_thread_id;
	ncol = graph->n;
	nrow = bandwidth + 1;
	*bandmatrix = Calloc(ncol * nrow, double);

#pragma omp parallel for private(i)
	for (i = 0; i < graph->n; i++) {
		int node = remap[i];
		int j;
		double val;

		GMRFLib_thread_id = id;

		val = Qfunc(i, i, Qfunc_arg);
		GMRFLib_STOP_IF_NAN_OR_INF(val, i, i);
		(*bandmatrix)[BIDX(0, node)] = val;

		for (j = 0; j < graph->nnbs[i]; j++) {
			int jj = graph->nbs[i][j];
			int nnode = remap[jj];

			if (nnode > node) {
				val = Qfunc(i, jj, Qfunc_arg);
				GMRFLib_STOP_IF_NAN_OR_INF(val, i, jj);
				(*bandmatrix)[BIDX(nnode - node, node)] = val;
			}
		}
	}

	if (GMRFLib_catch_error_for_inla) {
		if (nan_error) {
			return !GMRFLib_SUCCESS;
		}
	}

	if (0) {
		static int count = 0;
		char *fnm = NULL;

		GMRFLib_sprintf(&fnm, "Q-band-%1d-%1d.dat", count++, omp_get_thread_num());
		FILE *fp = fopen(fnm, "w");
		Free(fnm);
		assert(fp);
		FIXME("write Q-file");
		for (i = 0; i < graph->n; i++) {
			int node = remap[i];
			int j;

			fprintf(fp, "%d %d %.20f\n", i, i, (*bandmatrix)[BIDX(0, node)]);
			for (j = 0; j < graph->nnbs[i]; j++) {
				int jj = graph->nbs[i][j];
				int nnode = remap[jj];
				if (nnode > node) {
					fprintf(fp, "%d %d %.20f\n", i, jj, (*bandmatrix)[BIDX(nnode - node, node)]);
				}
			}
		}
		fclose(fp);
	}

	GMRFLib_thread_id = id;

	return GMRFLib_SUCCESS;
#undef BIDX
}

int GMRFLib_factorise_sparse_matrix_BAND(double *band, GMRFLib_fact_info_tp * finfo, GMRFLib_graph_tp * graph, int bandwidth)
{
	/*
	 * compute the factorisation of 'band' and overwrite it with the Cholesky-factorisation 
	 */

	int error = 0, nband, ldim, i, k;

	nband = bandwidth;
	ldim = bandwidth + 1;

	switch (GMRFLib_blas_level) {
	case BLAS_LEVEL2:
		dpbtf2_("L", &(graph->n), &nband, band, &ldim, &error, 1);
		break;
	case BLAS_LEVEL3:
		dpbtrf_("L", &(graph->n), &nband, band, &ldim, &error, 1);
		break;
	}
	if (error) {
		if (GMRFLib_catch_error_for_inla) {
			fprintf(stdout, "\n\t%s\n\tFunction: %s(), Line: %1d, Thread: %1d\n\tFail to factorize Q. I will try to fix it...\n\n",
				RCSId, __GMRFLib_FuncName, __LINE__, omp_get_thread_num());
			return GMRFLib_EPOSDEF;
		} else {
			GMRFLib_ERROR(GMRFLib_EPOSDEF);
		}
	}

	/*
	 * provide some info about the factorization 
	 */
	for (i = 0, k = 0; i < graph->n; i++)
		k += graph->nnbs[i];

	finfo->n = graph->n;				       /* size of Q */
	finfo->nnzero = k + graph->n;			       /* # non-zeros in Q */

	/*
	 * this is correct, as it will 'compute' 0's! 
	 */
	for (i = 0, k = 0; i < graph->n; i++)
		k += IMIN(graph->n, i + nband + 1) - (i + 1);

	finfo->nfillin = k - (finfo->nnzero - graph->n) / 2;   /* fillin in L */

	return GMRFLib_SUCCESS;
}

int GMRFLib_free_fact_sparse_matrix_BAND(double *bchol)
{
	Free(bchol);
	return GMRFLib_SUCCESS;
}

int GMRFLib_solve_lt_sparse_matrix_BAND(double *rhs, double *bchol, GMRFLib_graph_tp * graph, int *remap, int bandwidth)
{
	/*
	 * rhs in real world, bchol in mapped word
	 * 
	 * solve L^Tx=rhs, rhs is overwritten by the solution 
	 */
	int nband, ldim, stride = 1;

	nband = bandwidth;
	ldim = nband + 1;

	GMRFLib_convert_to_mapped(rhs, NULL, graph, remap);
	dtbsv_("L", "T", "N", &(graph->n), &nband, bchol, &ldim, rhs, &stride, 1, 1, 1);

	GMRFLib_convert_from_mapped(rhs, NULL, graph, remap);
	return GMRFLib_SUCCESS;
}

int GMRFLib_solve_llt_sparse_matrix_BAND(double *rhs, double *bchol, GMRFLib_graph_tp * graph, int *remap, int bandwidth)
{
	/*
	 * rhs in real world, bchol in mapped word
	 * 
	 * solve Q x=rhs, where Q=L L^T 
	 */
	int nband, ldim, stride = 1;

	nband = bandwidth;
	ldim = nband + 1;

	GMRFLib_convert_to_mapped(rhs, NULL, graph, remap);
	dtbsv_("L", "N", "N", &(graph->n), &nband, bchol, &ldim, rhs, &stride, 1, 1, 1);
	dtbsv_("L", "T", "N", &(graph->n), &nband, bchol, &ldim, rhs, &stride, 1, 1, 1);
	GMRFLib_convert_from_mapped(rhs, NULL, graph, remap);

	return GMRFLib_SUCCESS;
}

int GMRFLib_solve_l_sparse_matrix_BAND(double *rhs, double *bchol, GMRFLib_graph_tp * graph, int *remap, int bandwidth)
{
	/*
	 * rhs in real world, bchol in mapped word
	 * 
	 * solve Q x=rhs, where Q=L L^T 
	 */
	int nband, ldim, stride = 1;

	nband = bandwidth;
	ldim = nband + 1;

	GMRFLib_convert_to_mapped(rhs, NULL, graph, remap);
	dtbsv_("L", "N", "N", &(graph->n), &nband, bchol, &ldim, rhs, &stride, 1, 1, 1);
	GMRFLib_convert_from_mapped(rhs, NULL, graph, remap);

	return GMRFLib_SUCCESS;
}

int GMRFLib_solve_lt_sparse_matrix_special_BAND(double *rhs, double *bchol, GMRFLib_graph_tp * graph, int *remap, int bandwidth,
						int findx, int toindx, int remapped)
{
	/*
	 * rhs in real world, bchol in mapped world.  solve L^Tx=b backward only from rhs[findx] up to rhs[toindx].  note that
	 * findx and toindx is in mapped world.  if remapped, do not remap/remap-back the rhs before solving.
	 * 
	 */
	int nband, ldim, stride = 1, from, to;

	from = findx + 1;				       /* convert to fortran indxing */
	to = toindx + 1;				       /* convert to fortran indxing */
	nband = bandwidth;
	ldim = nband + 1;

	if (!remapped) {
		GMRFLib_convert_to_mapped(rhs, NULL, graph, remap);
	}
	dtbsvspecial_("L", "T", "N", &(graph->n), &nband, bchol, &ldim, rhs, &stride, &from, &to, 1, 1, 1);
	if (!remapped) {
		GMRFLib_convert_from_mapped(rhs, NULL, graph, remap);
	}

	return GMRFLib_SUCCESS;
}

int GMRFLib_solve_l_sparse_matrix_special_BAND(double *rhs, double *bchol, GMRFLib_graph_tp * graph, int *remap, int bandwidth, int findx,
					       int toindx, int remapped)
{
	/*
	 * rhs in real world, bchol in mapped world.  solve Lx=b backward only from rhs[findx] up to rhs[toindx].  note that
	 * findx and toindx is in mapped world.  if remapped, do not remap/remap-back the rhs before solving.
	 * 
	 */
	int nband, ldim, stride = 1, from, to;

	from = findx + 1;				       /* convert to fortran indxing */
	to = toindx + 1;				       /* convert to fortran indxing */
	nband = bandwidth;
	ldim = nband + 1;

	if (!remapped) {
		GMRFLib_convert_to_mapped(rhs, NULL, graph, remap);
	}
	dtbsvspecial_("L", "N", "N", &(graph->n), &nband, bchol, &ldim, rhs, &stride, &from, &to, 1, 1, 1);
	if (!remapped) {
		GMRFLib_convert_from_mapped(rhs, NULL, graph, remap);
	}

	return GMRFLib_SUCCESS;
}

int GMRFLib_solve_llt_sparse_matrix_special_BAND(double *rhs, double *bchol, GMRFLib_graph_tp * graph, int *remap, int bandwidth, int idx)
{
	/*
	 * rhs in real world, bchol in mapped word
	 * 
	 * solve Q x=rhs, where Q=L L^T 
	 */
	GMRFLib_ASSERT(rhs[idx] == 1.0, GMRFLib_ESNH);

	int nband, ldim, stride = 1, idxnew, from, to;

	nband = bandwidth;
	ldim = nband + 1;
	idxnew = remap[idx];
	rhs[idx] = 0.0;
	rhs[idxnew] = 1.0;
	from = idxnew + 1;
	to = graph->n;

	dtbsvspecial_("L", "N", "N", &(graph->n), &nband, bchol, &ldim, rhs, &stride, &from, &to, 1, 1, 1);
	dtbsv_("L", "T", "N", &(graph->n), &nband, bchol, &ldim, rhs, &stride, 1, 1, 1);
	GMRFLib_convert_from_mapped(rhs, NULL, graph, remap);

	if (0) {
		double *rrhs = Calloc(graph->n, double);
		rrhs[idx] = 1.0;
		GMRFLib_solve_llt_sparse_matrix_BAND(rrhs, bchol, graph, remap, bandwidth);
		for (int i = 0; i < graph->n; i++)
			fprintf(stderr, "%d %g %g %g\n", i, rhs[i], rrhs[i], rhs[i] - rrhs[i]);
	}

	return GMRFLib_SUCCESS;
}


int GMRFLib_comp_cond_meansd_BAND(double *cmean, double *csd, int indx, double *x, int remapped, double *bchol, GMRFLib_graph_tp * graph,
				  int *remap, int bandwidth)
{
	/*
	 * compute the conditonal mean and stdev for x[indx]|x[indx+1]...x[n-1] for the current value of x. if `remapped', then 
	 * x is assumed to be remapped for possible (huge) speedup when used repeately.  note: indx is still the user world!!!
	 * 
	 * example: approach 1,2,3 are equivalent (old style!)
	 * 
	 * (*GMRFLib_uniform_init)(seed); set_stdgauss(x); memcpy(z, x, graph->n*sizeof(double));
	 * 
	 * a1: gmrf_g_solve(x, bchol, graph);
	 * 
	 * a2: for(i=graph->n-1;i>=0;i--){ gmrf_g_compute_cmean_csd(&cmean, &csd, i, y, 0, bchol, graph); y[i] = cmean +
	 * csd*z[i];}
	 * 
	 * a3: for(i=graph->n-1;i>=0;i--){ gmrf_g_compute_cmean_csd(&cmean, &csd, i, y, 1, bchol, graph); y[graph->remap[i]] =
	 * cmean + csd*z[i];} 
	 */

	int nband, ldim, ii;

	nband = bandwidth;
	ldim = nband + 1;
	ii = remap[indx] + 1;				       /* fortran indxing */

	if (remapped)
		cmsd_(cmean, csd, &ii, &(graph->n), &nband, bchol, &ldim, x);
	else {
		GMRFLib_convert_to_mapped(x, NULL, graph, remap);
		cmsd_(cmean, csd, &ii, &(graph->n), &nband, bchol, &ldim, x);
		GMRFLib_convert_from_mapped(x, NULL, graph, remap);
	}
	return GMRFLib_SUCCESS;
}

int GMRFLib_log_determinant_BAND(double *logdet, double *bchol, GMRFLib_graph_tp * graph, int bandwidth)
{
	int ldim = bandwidth + 1, i;

	for (i = 0, *logdet = 0.0; i < graph->n; i++)
		*logdet += log(bchol[i * ldim]);
	*logdet *= 2;

	return GMRFLib_SUCCESS;
}

int GMRFLib_compute_Qinv_BAND(GMRFLib_problem_tp * problem, int storage)
{
	/*
	 * large parts of this code is copied from the TAUCS version. be aware.... 
	 */

#define BIDX(i,j)        ((j)+(i)*ldim)
#define LIDX(i,j)        BIDX(j, (i)-(j))
#define L(i,j)           Lmatrix[LIDX(i, j)]
#define Cov(i,j)         cov[LIDX(IMAX(i, j), IMIN(i, j))]
#define Cov_maxmin(i,j)  cov[LIDX(i, j)]		       /* save the MAX and the MIN */

	int i, j, k, kk, iii, jjj, bw, ldim, n, *inv_remap = NULL, *rremove = NULL, nrremove;
	double tmp, Lii_inv, value, *Lmatrix, *cov;

	map_id **Qinv_L = NULL;
	map_ii *mapping = NULL;

	bw = problem->sub_sm_fact.bandwidth;
	ldim = bw + 1;
	n = problem->sub_graph->n;
	Lmatrix = problem->sub_sm_fact.bchol;
	cov = Calloc(n * ldim, double);

	/*
	 * this block is copied from smtp-taucs.c
	 * 
	 * setup the hash-table for storing Qinv_L 
	 */
	Qinv_L = Calloc(n, map_id *);
#pragma omp parallel for private(i)
	for (i = 0; i < n; i++) {
		Qinv_L[i] = Calloc(1, map_id);
		map_id_init_hint(Qinv_L[i], ldim);
	}

	/*
	 * do the recursions 
	 */
	if (0) {
		/*
		 * plain and safe version, DO NOT DELETE!!!! 
		 */
		for (i = n - 1; i >= 0; i--) {
			Lii_inv = 1.0 / L(i, i);
			for (j = IMIN(i + bw, n - 1); j >= i; j--) {
				for (k = i + 1, tmp = 0.0; k <= IMIN(i + bw, n - 1); k++) {
					tmp += Cov(k, j) * L(k, i);
				}
				if (i == j) {
					Cov(i, j) = value = Lii_inv * (Lii_inv - tmp);
				} else {
					Cov(i, j) = value = -Lii_inv * tmp;
				}
				map_id_set(Qinv_L[i], j, value);
			}
		}
	} else {
		/*
		 * fast version. to understand this version: split the loop into two; k >= j and k < j, and also for i=j and
		 * i!=j. then expand the macros in the plain version.
		 * 
		 * TODO: it is possible to avoid storing all computed covariances twice, both in the hash-array and then also
		 * in 'cov'. there is something to be learned from the taucs-code of this function, where that part of Sigma
		 * that is needed is extracted from Qinv_L. However, this function is not ment to be used for large problems,
		 * so.... 
		 */
		for (i = n - 1; i >= 0; i--) {
			double *cov_offset, *Lmatrix_offset, *cov_set_offset;

			Lmatrix_offset = &Lmatrix[i * ldim - i];
			Lii_inv = 1.0 / L(i, i);

			cov_set_offset = &cov[i * bw];

			for (j = IMIN(i + bw, n - 1); j > i; j--) {
				tmp = 0.0;
				cov_offset = &cov[j * bw];
				for (k = IMAX(i + 1, j); k <= IMIN(i + bw, n - 1); k++) {
					tmp += cov_offset[k] * Lmatrix_offset[k];
				}
				cov_offset = &cov[j + (i + 1) * bw];
				for (k = i + 1, kk = 0; k < IMAX(i + 1, j); k++, kk += bw) {
					tmp += cov_offset[kk] * Lmatrix_offset[k];
				}
				cov_set_offset[j] = value = -Lii_inv * tmp;
				map_id_set(Qinv_L[i], j, value);
			}

			/*
			 * this is for j=i 
			 */
			{
				tmp = 0.0;

				cov_offset = &cov[i * bw];
				for (k = i + 1; k <= IMIN(i + bw, n - 1); k++) {
					tmp += cov_offset[k] * Lmatrix_offset[k];
				}
				cov_set_offset[i] = value = Lii_inv * (Lii_inv - tmp);
				map_id_set(Qinv_L[i], i, value);
			}
		}
	}

	/*
	 * ----------->>>> the rest of this routine is copied from smtp-taucs.c <<<<<--------------- 
	 */

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
	Free(cov);
	Free(inv_remap);
	Free(rremove);

#undef BIDX
#undef LIDX
#undef L
#undef Cov
	return GMRFLib_SUCCESS;
}

/* 
   from here on: undocumented features
*/

/* 
   from here is for internal use only. not documented
*/
int GMRFLib_bitmap_factorisation_BAND__intern(const char *filename, double *band, GMRFLib_graph_tp * graph, int *remap, int bandwidth)
{
#define BIDX(i,j) ((i)+(j)*ldim)			       /* band index'ing */

#define NBitsInByte 8
#define SETBIT(im, jm, m, N) {						\
		int local_im = (int)((im) * reduce_factor);		\
		int local_jm = (int)((jm) * reduce_factor);		\
		if (GMRFLib_bitmap_swap){				\
			int itmp = local_im;				\
			local_im = N-1-local_jm;			\
			local_jm = itmp;				\
		}							\
		int ii = (local_im)/NBitsInByte;			\
		GMRFLib_setbit(&bitmap[ ii+(local_jm)*(m)], (unsigned int) (NBitsInByte-1-((local_im)-ii*NBitsInByte))); \
	}

	int i, j, n = graph->n, N, m;
	double reduce_factor;
	unsigned char *bitmap;
	FILE *fp;

	int nband = bandwidth;
	int ldim = bandwidth + 1;

	if (GMRFLib_bitmap_max_dimension > 0 && n > GMRFLib_bitmap_max_dimension) {
		N = GMRFLib_bitmap_max_dimension;
		reduce_factor = (double) N / (double) n;
	} else {
		N = n;
		reduce_factor = 1.0;
	}

	m = N / NBitsInByte;
	if (m * NBitsInByte != N)
		m++;
	bitmap = Calloc(m * N, unsigned char);

	for (i = 0; i < graph->n; i++) {
		for (j = i; j < IMIN(i + nband + 1, graph->n); j++) {
			if (!ISZERO(band[BIDX(j - i, i)])) {
				SETBIT(i, j, m, N);
			}
		}
	}

	fp = fopen(filename, "w");
	if (fp) {
		fprintf(fp, "P4\n%1d %1d\n", N, N);
		for (i = 0; i < N; i++) {
			fwrite(&bitmap[i * m], (unsigned int) m, 1, fp);
		}
		fclose(fp);
	} else {
		GMRFLib_ERROR(GMRFLib_EOPENFILE);
	}
	Free(bitmap);
	return GMRFLib_SUCCESS;
#undef SETBIT
#undef NBitsInByte


#undef SETBIT
#undef NBitsInByte
#undef BIDX
}

int GMRFLib_bitmap_factorisation_BAND(const char *filename_body, double *band, GMRFLib_graph_tp * graph, int *remap, int bandwidth)
{
	/*
	 * create a bitmap-file of the factorization 
	 */
	char *filename;

	GMRFLib_EWRAP0(GMRFLib_sprintf(&filename, "%s_L.pbm", (filename_body ? filename_body : "band_L")));
	GMRFLib_EWRAP0(GMRFLib_bitmap_factorisation_BAND__intern(filename, band, graph, remap, bandwidth));
	Free(filename);

	return GMRFLib_SUCCESS;
}

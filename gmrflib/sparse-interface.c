
/* GMRFLib-sparse-interface.c
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
  \file sparse-interface.c
  \brief Unified interface to the sparse-matrix libraries
*/

#include "GMRFLib/GMRFLib.h"
#include "GMRFLib/GMRFLibP.h"

/*!
  \brief Compute the reordering
*/
int GMRFLib_compute_reordering(GMRFLib_sm_fact_tp *sm_fact, GMRFLib_graph_tp *graph, GMRFLib_global_node_tp *gn)
{
	GMRFLib_ENTER_ROUTINE;
	GMRFLib_global_node_tp lgn, *gn_ptr = NULL;

	if (gn) {
		/*
		 * then this defines the global node definition 
		 */
		gn_ptr = (GMRFLib_global_node_tp *) gn;
	} else {
		lgn = GMRFLib_global_node;
		gn_ptr = &lgn;
	}

	switch (GMRFLib_reorder) {
	case GMRFLib_REORDER_DEFAULT:
	{
		/*
		 * this choice depends on the sparse-solver 
		 */
		switch (sm_fact->smtp) {
		case GMRFLib_SMTP_BAND:
		{
			GMRFLib_EWRAP1(GMRFLib_compute_reordering_BAND(&(sm_fact->remap), graph));
		}
			break;

		case GMRFLib_SMTP_TAUCS:
		{
			GMRFLib_EWRAP1(GMRFLib_compute_reordering_TAUCS(&(sm_fact->remap), graph, GMRFLib_reorder, gn_ptr));
		}
			break;

		case GMRFLib_SMTP_PARDISO:
		{
			if (sm_fact->PARDISO_fact == NULL) {
				GMRFLib_pardiso_init(&(sm_fact->PARDISO_fact));
			}
			GMRFLib_pardiso_reorder(sm_fact->PARDISO_fact, graph);
			sm_fact->remap = Calloc(graph->n, int);
			Memcpy((void *) sm_fact->remap, (void *) sm_fact->PARDISO_fact->pstore[GMRFLib_PSTORE_TNUM_REF]->perm,
			       graph->n * sizeof(int));
		}
			break;

		default:
			GMRFLib_ASSERT(1 == 0, GMRFLib_ESNH);
			break;
		}
	}
		break;

	case GMRFLib_REORDER_PARDISO:
	{
		switch (sm_fact->smtp) {
		case GMRFLib_SMTP_PARDISO:		       /* same code as above */
		{
			if (sm_fact->PARDISO_fact == NULL) {
				GMRFLib_pardiso_init(&(sm_fact->PARDISO_fact));
			}
			GMRFLib_pardiso_reorder(sm_fact->PARDISO_fact, graph);
			sm_fact->remap = Calloc(graph->n, int);
			Memcpy((void *) sm_fact->remap, (void *) sm_fact->PARDISO_fact->pstore[GMRFLib_PSTORE_TNUM_REF]->perm,
			       graph->n * sizeof(int));
		}
			break;

		default:
			GMRFLib_ASSERT(1 == 0, GMRFLib_ESNH);
			break;
		}
	}
		break;

	case GMRFLib_REORDER_BAND:
		GMRFLib_EWRAP1(GMRFLib_compute_reordering_BAND(&(sm_fact->remap), graph));
		break;

		/*
		 * all the remaining ones are treated by the _TAUCS routine 
		 */
	case GMRFLib_REORDER_IDENTITY:
	case GMRFLib_REORDER_REVERSE_IDENTITY:
	case GMRFLib_REORDER_METIS:
	case GMRFLib_REORDER_GENMMD:
	case GMRFLib_REORDER_AMD:
	case GMRFLib_REORDER_AMDC:
	case GMRFLib_REORDER_AMDBAR:
	case GMRFLib_REORDER_AMDBARC:
	case GMRFLib_REORDER_MD:
	case GMRFLib_REORDER_MMD:
	{
		GMRFLib_EWRAP1(GMRFLib_compute_reordering_TAUCS(&(sm_fact->remap), graph, GMRFLib_reorder, gn_ptr));
	}
		break;
	default:
		P(GMRFLib_reorder);
		GMRFLib_ASSERT(0 == 1, GMRFLib_ESNH);
	}

	if (sm_fact->remap) {				       /* need this still for the wa-routines. FIXME */
		GMRFLib_EWRAP1(GMRFLib_graph_comp_bw(&(sm_fact->bandwidth), graph, sm_fact->remap));
	} else {
		sm_fact->bandwidth = -1;
	}

	GMRFLib_LEAVE_ROUTINE;

	return GMRFLib_SUCCESS;
}

/*!
  \brief Free the reordering
*/
int GMRFLib_free_reordering(GMRFLib_sm_fact_tp *sm_fact)
{
	if (sm_fact) {
		Free(sm_fact->remap);
		sm_fact->bandwidth = 0;
	}
	return GMRFLib_SUCCESS;
}

/*
  \brief Build a sparse matrix
*/
int GMRFLib_build_sparse_matrix(int thread_id, GMRFLib_sm_fact_tp *sm_fact, GMRFLib_Qfunc_tp *Qfunc, void *Qfunc_arg, GMRFLib_graph_tp *graph)
{
	GMRFLib_ENTER_ROUTINE;
	int ret;

	switch (sm_fact->smtp) {
	case GMRFLib_SMTP_BAND:
	{
		ret = GMRFLib_build_sparse_matrix_BAND(thread_id, &(sm_fact->bchol), Qfunc, Qfunc_arg, graph, sm_fact->remap, sm_fact->bandwidth);
		if (ret != GMRFLib_SUCCESS) {
			return ret;
		}
	}
		break;

	case GMRFLib_SMTP_TAUCS:
	{
		ret = GMRFLib_build_sparse_matrix_TAUCS(thread_id, &(sm_fact->TAUCS_L), Qfunc, Qfunc_arg, graph, sm_fact->remap);
		if (ret != GMRFLib_SUCCESS) {
			return ret;
		}
	}
		break;

	case GMRFLib_SMTP_PARDISO:
	{
		if (sm_fact->PARDISO_fact == NULL) {
			GMRFLib_pardiso_init(&(sm_fact->PARDISO_fact));
			GMRFLib_pardiso_reorder(sm_fact->PARDISO_fact, graph);
		}
		ret = GMRFLib_pardiso_build(thread_id, sm_fact->PARDISO_fact, graph, Qfunc, Qfunc_arg);
		if (ret != GMRFLib_SUCCESS) {
			return ret;
		}
	}
		break;

	default:
		GMRFLib_ASSERT(1 == 0, GMRFLib_ESNH);
		break;
	}

	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}

/*!
  \brief Factorise a sparse matrix
*/
int GMRFLib_factorise_sparse_matrix(GMRFLib_sm_fact_tp *sm_fact, GMRFLib_graph_tp *graph)
{
	int ret;
	GMRFLib_ENTER_ROUTINE;

	switch (sm_fact->smtp) {
	case GMRFLib_SMTP_BAND:
	{
		ret = GMRFLib_factorise_sparse_matrix_BAND(sm_fact->bchol, &(sm_fact->finfo), graph, sm_fact->bandwidth);
		if (ret != GMRFLib_SUCCESS) {
			return ret;
		}
	}
		break;

	case GMRFLib_SMTP_TAUCS:
	{
		ret = GMRFLib_factorise_sparse_matrix_TAUCS(&(sm_fact->TAUCS_L), &(sm_fact->TAUCS_symb_fact), &(sm_fact->finfo),
							    &(sm_fact->TAUCS_L_inv_diag));
		if (ret != GMRFLib_SUCCESS) {
			return ret;
		}
	}
		break;

	case GMRFLib_SMTP_PARDISO:
	{
		ret = GMRFLib_pardiso_chol(sm_fact->PARDISO_fact);
		if (ret != GMRFLib_SUCCESS) {
			return ret;
		}
	}
		break;

	default:
		GMRFLib_ASSERT(1 == 0, GMRFLib_ESNH);
		break;
	}

	GMRFLib_LEAVE_ROUTINE;

	return GMRFLib_SUCCESS;
}

/*!
  \brief Free a factorisation of a sparse matrix
*/
int GMRFLib_free_fact_sparse_matrix(GMRFLib_sm_fact_tp *sm_fact)
{
	if (sm_fact) {
		switch (sm_fact->smtp) {
		case GMRFLib_SMTP_BAND:
		{
			GMRFLib_EWRAP1(GMRFLib_free_fact_sparse_matrix_BAND(sm_fact->bchol));
			sm_fact->bchol = NULL;
		}
			break;

		case GMRFLib_SMTP_TAUCS:
		{
			GMRFLib_EWRAP1(GMRFLib_free_fact_sparse_matrix_TAUCS
				       (sm_fact->TAUCS_L, sm_fact->TAUCS_L_inv_diag, sm_fact->TAUCS_symb_fact));
			sm_fact->TAUCS_L = NULL;
			sm_fact->TAUCS_symb_fact = NULL;
		}
			break;

		case GMRFLib_SMTP_PARDISO:
		{
			if (sm_fact->PARDISO_fact) {
				GMRFLib_EWRAP1(GMRFLib_pardiso_free(&(sm_fact->PARDISO_fact)));
				sm_fact->PARDISO_fact = NULL;
			}
		}
			break;

		default:
			GMRFLib_ASSERT(1 == 0, GMRFLib_ESNH);
			break;
		}
	}

	return GMRFLib_SUCCESS;
}

/*!
  \brief Solve \f$Lx=b\f$
*/
int GMRFLib_solve_l_sparse_matrix(double *rhs, int nrhs, GMRFLib_sm_fact_tp *sm_fact, GMRFLib_graph_tp *graph)
{
	/*
	 * rhs in real world. solve L x=rhs, rhs is overwritten by the solution 
	 */
	GMRFLib_ENTER_ROUTINE;

	switch (sm_fact->smtp) {
	case GMRFLib_SMTP_BAND:
	{
		omp_set_num_threads(GMRFLib_openmp->max_threads_inner);
#pragma omp parallel for num_threads(GMRFLib_openmp->max_threads_inner)
		for (int i = 0; i < nrhs; i++) {
			GMRFLib_solve_l_sparse_matrix_BAND(&rhs[i * graph->n], sm_fact->bchol, graph, sm_fact->remap, sm_fact->bandwidth);
		}
	}
		break;

	case GMRFLib_SMTP_TAUCS:
	{
		omp_set_num_threads(GMRFLib_openmp->max_threads_inner);
#pragma omp parallel for num_threads(GMRFLib_openmp->max_threads_inner)
		for (int i = 0; i < nrhs; i++) {
			GMRFLib_solve_l_sparse_matrix_TAUCS(&rhs[i * graph->n], sm_fact->TAUCS_L, graph, sm_fact->remap);
		}
	}
		break;

	case GMRFLib_SMTP_PARDISO:
	{
		GMRFLib_pardiso_solve_L(sm_fact->PARDISO_fact, rhs, rhs, nrhs);
	}
		break;

	default:
		GMRFLib_ERROR(GMRFLib_ESNH);
		break;
	}

	GMRFLib_LEAVE_ROUTINE;

	return GMRFLib_SUCCESS;
}

/*!
  \brief Solve \f$L^Tx=b\f$
*/
int GMRFLib_solve_lt_sparse_matrix(double *rhs, int nrhs, GMRFLib_sm_fact_tp *sm_fact, GMRFLib_graph_tp *graph)
{
	/*
	 * rhs in real world. solve L^Tx=rhs, rhs is overwritten by the solution 
	 */
	GMRFLib_ENTER_ROUTINE;

	switch (sm_fact->smtp) {
	case GMRFLib_SMTP_BAND:
	{
		omp_set_num_threads(GMRFLib_openmp->max_threads_inner);
#pragma omp parallel for num_threads(GMRFLib_openmp->max_threads_inner)
		for (int i = 0; i < nrhs; i++) {
			GMRFLib_solve_lt_sparse_matrix_BAND(&rhs[i * graph->n], sm_fact->bchol, graph, sm_fact->remap, sm_fact->bandwidth);
		}
	}
		break;

	case GMRFLib_SMTP_TAUCS:
	{
		omp_set_num_threads(GMRFLib_openmp->max_threads_inner);
#pragma omp parallel for num_threads(GMRFLib_openmp->max_threads_inner)
		for (int i = 0; i < nrhs; i++) {
			GMRFLib_solve_lt_sparse_matrix_TAUCS(&rhs[i * graph->n], sm_fact->TAUCS_L, graph, sm_fact->remap);
		}
	}
		break;

	case GMRFLib_SMTP_PARDISO:
	{
		GMRFLib_pardiso_solve_LT(sm_fact->PARDISO_fact, rhs, rhs, nrhs);
	}
		break;

	default:
		GMRFLib_ERROR(GMRFLib_ESNH);
		break;
	}

	GMRFLib_LEAVE_ROUTINE;

	return GMRFLib_SUCCESS;
}

/*!
  \brief Solve \f$LL^Tx=b\f$  or \f$Qx=b\f$
*/
int GMRFLib_solve_llt_sparse_matrix(double *rhs, int nrhs, GMRFLib_sm_fact_tp *sm_fact, GMRFLib_graph_tp *graph)
{
	/*
	 * rhs in real world. solve Q x=rhs, where Q=L L^T 
	 */
	GMRFLib_ENTER_ROUTINE;

	if (sm_fact->smtp == GMRFLib_SMTP_BAND) {
		omp_set_num_threads(GMRFLib_openmp->max_threads_inner);
#pragma omp parallel for num_threads(GMRFLib_openmp->max_threads_inner)
		for (int i = 0; i < nrhs; i++) {
			GMRFLib_solve_llt_sparse_matrix_BAND(&rhs[i * graph->n], sm_fact->bchol, graph, sm_fact->remap, sm_fact->bandwidth);
		}
	} else if (sm_fact->smtp == GMRFLib_SMTP_TAUCS) {
		int ntt = -1;
		if (omp_get_level() == 0) {
			ntt = GMRFLib_PARDISO_MAX_NUM_THREADS();
		} else {
			ntt = GMRFLib_openmp->max_threads_inner;
		}

		int numt_save = omp_get_max_threads();
		int reset_num_threads = 0;

		if (nrhs <= ntt * 4) {
			if (nrhs > 1) {
				omp_set_num_threads(IMIN(nrhs, ntt));
				reset_num_threads = 1;
#pragma omp parallel for
				for (int i = 0; i < nrhs; i++) {
					GMRFLib_solve_llt_sparse_matrix_TAUCS(&rhs[i * graph->n], sm_fact->TAUCS_L, graph, sm_fact->remap);
				}
			} else {
				GMRFLib_solve_llt_sparse_matrix_TAUCS(rhs, sm_fact->TAUCS_L, graph, sm_fact->remap);
			}
		} else {
			// much of the same code as in the smtp-pardiso.c and solve_core function
			int nt = 1;
			int nsolve;
			int nblock;
			int block_nrhs;
			div_t d;
			int S_nrhs_max = GMRFLib_pardiso_get_nrhs();

			if (nrhs > 1) {
				if (GMRFLib_openmp->adaptive && omp_get_level() == 0) {
					// this is the exception of the rule, as we want to run this in parallel if we are in adaptive model and
					// level=0.
					nt = GMRFLib_PARDISO_MAX_NUM_THREADS();
				} else {
					nt = GMRFLib_openmp->max_threads_inner;
				}
			}
			if (nrhs == 1) {
				// in this case we always set block_nrhs=1 and nt=1
				block_nrhs = 1;
				nt = 1;
			} else if (S_nrhs_max < 0) {
				// this is the adaptive choice
				d = div(nrhs, nt);
				if (nrhs >= nt) {
					// if this is true, we need to divide the work between the nt cores
					block_nrhs = d.quot + (d.rem != 0);
				} else {
					// else we use one core pr rhs
					block_nrhs = 1;
					nt = IMIN(nt, nrhs);
				}
			} else if (nrhs <= S_nrhs_max) {
				// then we do all of them in one block
				block_nrhs = nrhs;
				nt = 1;
			} else {
				// S_nrhs_max define the max nrhs, so then we divide the work
				block_nrhs = S_nrhs_max;
				d = div(nrhs, block_nrhs);
				nt = IMIN(nt, d.quot + (d.rem != 0));
			}

			if (nt > 1) {
				omp_set_num_threads(nt);
				reset_num_threads = 1;
			}

			d = div(nrhs, block_nrhs);
			nblock = d.quot;
			nsolve = nblock + (d.rem != 0);

			if (0) {
				printf("nblock %d nsolve %d nrhs %d S_nrhs_max %d block_nrhs %d nt %d\n",
				       nblock, nsolve, nrhs, S_nrhs_max, block_nrhs, nt);
			}
#pragma omp parallel for num_threads(nt) if (nt > 1)
			for (int k = 0; k < nsolve; k++) {
				int offset = k * graph->n * block_nrhs;
				int local_nrhs = (k < nblock ? block_nrhs : (int) d.rem);
				GMRFLib_solve_llt_sparse_matrix2_TAUCS(rhs + offset, sm_fact->TAUCS_L, graph, sm_fact->remap, local_nrhs);
			}
		}

		if (reset_num_threads) {
			omp_set_num_threads(numt_save);
		}
	} else if (sm_fact->smtp == GMRFLib_SMTP_PARDISO) {
		GMRFLib_EWRAP1(GMRFLib_pardiso_solve_LLT(sm_fact->PARDISO_fact, rhs, rhs, nrhs));
	} else {
		GMRFLib_ERROR(GMRFLib_ESNH);
	}

	GMRFLib_LEAVE_ROUTINE;

	return GMRFLib_SUCCESS;
}

int GMRFLib_solve_llt_sparse_matrix_special(double *rhs, GMRFLib_sm_fact_tp *sm_fact, GMRFLib_graph_tp *graph, int idx)
{
	/*
	 * rhs in real world. solve Q x=rhs, where Q=L L^T. BUT, here we know that rhs is 0 execpt for a 1 at index idx.
	 */
	GMRFLib_ENTER_ROUTINE;

	switch (sm_fact->smtp) {
	case GMRFLib_SMTP_BAND:
	{
		GMRFLib_EWRAP1(GMRFLib_solve_llt_sparse_matrix_special_BAND(rhs, sm_fact->bchol, graph, sm_fact->remap, sm_fact->bandwidth, idx));
	}
		break;

	case GMRFLib_SMTP_TAUCS:
	{
		GMRFLib_EWRAP1(GMRFLib_solve_llt_sparse_matrix_special_TAUCS
			       (rhs, sm_fact->TAUCS_L, sm_fact->TAUCS_L_inv_diag, graph, sm_fact->remap, idx));
	}
		break;

	case GMRFLib_SMTP_PARDISO:
	{
		GMRFLib_pardiso_solve_LLT(sm_fact->PARDISO_fact, rhs, rhs, 1);
	}
		break;

	default:
		GMRFLib_ERROR(GMRFLib_ESNH);
		break;
	}

	GMRFLib_LEAVE_ROUTINE;

	return GMRFLib_SUCCESS;
}

/*!
  \brief Solve \f$L^Tx=b\f$ for indices in an interval
*/
int GMRFLib_solve_lt_sparse_matrix_special(double *rhs, GMRFLib_sm_fact_tp *sm_fact, GMRFLib_graph_tp *graph, int findx, int toindx, int remapped)
{
	/*
	 * rhs in real world, bchol in mapped world. solve L^Tx=b backward only from rhs[findx] up to rhs[toindx]. note that
	 * findx and toindx is in mapped world. if remapped, do not remap/remap-back the rhs before solving.
	 */
	GMRFLib_ENTER_ROUTINE;
	switch (sm_fact->smtp) {
	case GMRFLib_SMTP_BAND:
	{
		GMRFLib_EWRAP0(GMRFLib_solve_lt_sparse_matrix_special_BAND
			       (rhs, sm_fact->bchol, graph, sm_fact->remap, sm_fact->bandwidth, findx, toindx, remapped));
	}
		break;

	case GMRFLib_SMTP_TAUCS:
	{
		GMRFLib_EWRAP0(GMRFLib_solve_lt_sparse_matrix_special_TAUCS(rhs, sm_fact->TAUCS_L, graph, sm_fact->remap, findx, toindx, remapped));
	}
		break;

	case GMRFLib_SMTP_PARDISO:
	{
		GMRFLib_pardiso_solve_LT(sm_fact->PARDISO_fact, rhs, rhs, 1);
	}
		break;

	default:
	{
		GMRFLib_ERROR(GMRFLib_ESNH);
	}
		break;
	}
	GMRFLib_LEAVE_ROUTINE;

	return GMRFLib_SUCCESS;
}

/*!
  \brief Solve \f$Lx=b\f$ for indices in an interval
*/

int GMRFLib_solve_l_sparse_matrix_special(double *rhs, GMRFLib_sm_fact_tp *sm_fact, GMRFLib_graph_tp *graph, int findx, int toindx, int remapped)
{
	/*
	 * rhs in real world, bchol in mapped world. solve Lx=b backward only from rhs[findx] up to rhs[toindx]. note that
	 * findx and toindx is in mapped world. if remapped, do not remap/remap-back the rhs before solving.
	 */
	GMRFLib_ENTER_ROUTINE;
	switch (sm_fact->smtp) {
	case GMRFLib_SMTP_BAND:
	{
		GMRFLib_EWRAP0(GMRFLib_solve_l_sparse_matrix_special_BAND
			       (rhs, sm_fact->bchol, graph, sm_fact->remap, sm_fact->bandwidth, findx, toindx, remapped));
	}
		break;

	case GMRFLib_SMTP_TAUCS:
	{
		GMRFLib_EWRAP0(GMRFLib_solve_l_sparse_matrix_special_TAUCS(rhs, sm_fact->TAUCS_L, graph, sm_fact->remap, findx, toindx, remapped));
	}
		break;

	case GMRFLib_SMTP_PARDISO:
	{
		if (remapped) {
			GMRFLib_pardiso_perm(rhs, 1, sm_fact->PARDISO_fact);
		}
		GMRFLib_pardiso_solve_L(sm_fact->PARDISO_fact, rhs, rhs, 1);
	}
		break;

	default:
		GMRFLib_ERROR(GMRFLib_ESNH);
		break;
	}

	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}

/*!
  \brief Compute the log determininant of \f$Q\f$
*/
int GMRFLib_log_determinant(double *logdet, GMRFLib_sm_fact_tp *sm_fact, GMRFLib_graph_tp *graph)
{
	switch (sm_fact->smtp) {
	case GMRFLib_SMTP_BAND:
	{
		GMRFLib_EWRAP0(GMRFLib_log_determinant_BAND(logdet, sm_fact->bchol, graph, sm_fact->bandwidth));
	}
		break;

	case GMRFLib_SMTP_TAUCS:
	{
		GMRFLib_EWRAP0(GMRFLib_log_determinant_TAUCS(logdet, sm_fact->TAUCS_L));
	}
		break;

	case GMRFLib_SMTP_PARDISO:
	{
		*logdet = GMRFLib_pardiso_logdet(sm_fact->PARDISO_fact);
	}
		break;

	default:
		GMRFLib_ERROR(GMRFLib_ESNH);
		break;
	}

	return GMRFLib_SUCCESS;
}

/*!
  \brief Compute conditional mean and standard deviation of \f$x[i]\f$ conditioned on {\f$x[j]\f$}
    for \f$j>i\f$
*/
int GMRFLib_comp_cond_meansd(double *cmean, double *csd, int indx, double *x, int remapped, GMRFLib_sm_fact_tp *sm_fact, GMRFLib_graph_tp *graph)
{
	GMRFLib_ENTER_ROUTINE;
	switch (sm_fact->smtp) {
	case GMRFLib_SMTP_BAND:
	{
		GMRFLib_EWRAP1(GMRFLib_comp_cond_meansd_BAND
			       (cmean, csd, indx, x, remapped, sm_fact->bchol, graph, sm_fact->remap, sm_fact->bandwidth));
	}
		break;

	case GMRFLib_SMTP_TAUCS:
	{
		GMRFLib_EWRAP1(GMRFLib_comp_cond_meansd_TAUCS(cmean, csd, indx, x, remapped, sm_fact->TAUCS_L, graph, sm_fact->remap));
	}
		break;

	case GMRFLib_SMTP_PARDISO:
	{
		assert(0 == 1);
	}
		break;

	default:
		GMRFLib_ERROR(GMRFLib_ESNH);
		break;
	}
	GMRFLib_LEAVE_ROUTINE;

	return GMRFLib_SUCCESS;
}

/*!
  \brief Produce a bitmap of the Cholesky triangle in the portable bitmap (pbm) format
*/
int GMRFLib_bitmap_factorisation(const char *filename_body, GMRFLib_sm_fact_tp *sm_fact, GMRFLib_graph_tp *graph)
{
	switch (sm_fact->smtp) {
	case GMRFLib_SMTP_BAND:
	{
		GMRFLib_EWRAP1(GMRFLib_bitmap_factorisation_BAND(filename_body, sm_fact->bchol, graph, sm_fact->remap, sm_fact->bandwidth));
	}
		break;

	case GMRFLib_SMTP_TAUCS:
	{
		GMRFLib_EWRAP1(GMRFLib_bitmap_factorisation_TAUCS(filename_body, sm_fact->TAUCS_L));
	}
		break;

	case GMRFLib_SMTP_PARDISO:
	{
		GMRFLib_EWRAP1(GMRFLib_pardiso_bitmap());
	}
		break;

	default:
		GMRFLib_ERROR(GMRFLib_ESNH);
		break;
	}
	return GMRFLib_SUCCESS;
}

/*!
  \brief Wrapper for computing the (structural) inverse of \c Q.
*/
int GMRFLib_compute_Qinv(void *problem)
{
	GMRFLib_problem_tp *p = (GMRFLib_problem_tp *) problem;
	GMRFLib_ENTER_ROUTINE;

	switch (p->sub_sm_fact.smtp) {
	case GMRFLib_SMTP_BAND:
	{
		GMRFLib_EWRAP0(GMRFLib_compute_Qinv_BAND(p));
	}
		break;

	case GMRFLib_SMTP_TAUCS:
	{
		GMRFLib_EWRAP0(GMRFLib_compute_Qinv_TAUCS(p));
	}
		break;

	case GMRFLib_SMTP_PARDISO:
	{
		GMRFLib_pardiso_Qinv_INLA(p);
	}
		break;

	default:
		GMRFLib_ERROR(GMRFLib_ESNH);
		break;
	}
	GMRFLib_LEAVE_ROUTINE;
	return GMRFLib_SUCCESS;
}

/*!
  \brief Return \c GMRFLib_TRUE or \c GMRFLib_FALSE if smtp is valid
*/
int GMRFLib_valid_smtp(int smtp)
{
	if ((smtp == GMRFLib_SMTP_BAND) || (smtp == GMRFLib_SMTP_TAUCS)) {
		return GMRFLib_TRUE;
	} else {
		return GMRFLib_FALSE;
	}
}

/*! 
  \brief Return the name of a reordering
 */
const char *GMRFLib_reorder_name(GMRFLib_reorder_tp r)
{
	switch (r) {
	case GMRFLib_REORDER_DEFAULT:
		return "default";
	case GMRFLib_REORDER_IDENTITY:
		return "identity";
	case GMRFLib_REORDER_REVERSE_IDENTITY:
		return "reverseidentity";
	case GMRFLib_REORDER_BAND:
		return "band";
	case GMRFLib_REORDER_METIS:
		return "metis";
	case GMRFLib_REORDER_GENMMD:
		return "genmmd";
	case GMRFLib_REORDER_AMD:
		return "amd";
	case GMRFLib_REORDER_MD:
		return "md";
	case GMRFLib_REORDER_MMD:
		return "mmd";
	case GMRFLib_REORDER_AMDBAR:
		return "amdbar";
	case GMRFLib_REORDER_AMDC:
		return "amdc";
	case GMRFLib_REORDER_AMDBARC:
		return "amdbarc";
	case GMRFLib_REORDER_PARDISO:
		return "pardiso";
	default:
		fprintf(stderr, "\n\t*** ERROR *** Reordering [%d] not defined.\n", r);
		GMRFLib_ASSERT_RETVAL(0 == 1, GMRFLib_EPARAMETER, "(unknown reording)");
	}

	return "(unknown reording)";
}

/*! 
  \brief Return the id of the reordering from the name. 
 */
int GMRFLib_reorder_id(const char *name)
{
	if (!strcasecmp(name, "default") || !strcasecmp(name, "auto"))
		return GMRFLib_REORDER_DEFAULT;
	else if (!strcasecmp(name, "identity"))
		return GMRFLib_REORDER_IDENTITY;
	else if (!strcasecmp(name, "reverseidentity"))
		return GMRFLib_REORDER_REVERSE_IDENTITY;
	else if (!strcasecmp(name, "band"))
		return GMRFLib_REORDER_BAND;
	else if (!strcasecmp(name, "metis"))
		return GMRFLib_REORDER_METIS;
	else if (!strcasecmp(name, "genmmd"))
		return GMRFLib_REORDER_GENMMD;
	else if (!strcasecmp(name, "amd"))
		return GMRFLib_REORDER_AMD;
	else if (!strcasecmp(name, "amdbar"))
		return GMRFLib_REORDER_AMDBAR;
	else if (!strcasecmp(name, "md"))
		return GMRFLib_REORDER_MD;
	else if (!strcasecmp(name, "mmd"))
		return GMRFLib_REORDER_MMD;
	else if (!strcasecmp(name, "amdc"))
		return GMRFLib_REORDER_AMDC;
	else if (!strcasecmp(name, "amdbarc"))
		return GMRFLib_REORDER_AMDBARC;
	else if (!strcasecmp(name, "pardiso"))
		return GMRFLib_REORDER_PARDISO;
	else {
		fprintf(stderr, "\n\t*** ERROR *** Reordering [%s] not defined.\n", name);
		GMRFLib_ASSERT_RETVAL(0 == 1, GMRFLib_EPARAMETER, -1);
	}

	return -1;
}
